#####�J���e���p������̂��߂̐��ݕϐ����f�����O#####
library(MASS)
library(survival)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

set.seed(3181)

####�f�[�^�̓ǂݍ��݂ƃf�[�^���T���v�����O####
access_all <- read.csv("C:/Statistics/dat/access_log.csv", fileEncoding = "Shift-JIS")
company_all <- read.csv("C:/Statistics/dat/company_without_act_id.csv", fileEncoding = "Shift-JIS")
company_vcs <- read.csv("C:/Statistics/dat/company_with_act_id.csv", fileEncoding = "Shift-JIS")
access_all1 <- access_all

index.user1 <- subset(names(sort(table(access_all1$user_id), decreasing=TRUE)), 
                      sort(table(access_all1$user_id), decreasing=TRUE) < 500)
access_all <- access_all1[access_all1$user_id %in% index.user1, ]
sort(table(access_all$user_id), decreasing=TRUE)


#�s�v�ϐ�������
gv <- setdiff(colnames(access_all), c("u_regist_date", "url", "user_agent", "purpose1", "purpose2",
                                      "purpose3", "purpose4", "purpose5", "purpose6"))   #�s�v�ϐ�������


access_data <- access_all[, gv]


##���[�U�[id�������_���T���v�����O(5000�l)
n <- 10000
id_list <- unique(access_data$user_id)
id_sample <- sample(id_list, n)

#�����_���T���v�����O�����s
access_data_s <- access_data[access_data$user_id %in% id_sample, ]


####�f�[�^�̃N�����W���O####
####ID�̍쐬####
##user_id�𐔒l�ɕύX
td_id <- as.numeric(access_data_s$td_client_id)
user_id_o <- as.character(access_data_s$user_id)
unique_id <- unique(user_id_o)

#user_id��A�ԂɕύX
user_id <- array(0, dim=c(length(user_id_o)))
for(i in 1:length(unique_id)){
  print(i)
  user_id[user_id_o==unique_id[i]] <- i
}

##session_id�𐔒l�ɕύX
session_id_o <- as.numeric(access_data_s$session_id)
unique_id <- unique(user_id)

#session_id��A�ԂɕύX
session_id <- array(0, dim=c(length(session_id_o)))
for(i in 1:length(unique_id)){
  print(i)
  unique_session <- unique(session_id_o[user_id==unique_id[i]])
  
  #���[�U�[id���ƂɃZ�b�V����id�̘A�Ԃ�����
  for(j in 1:length(unique_session)){
    session_id[session_id_o==unique_session[j]] <- j
  }
}


#ID�̌����ƌ�������ID�̊m�F
ID <- data.frame(no=1:length(user_id), user_id, session_id)
ID_watch <- data.frame(user_id, session_id, user_id_o, session_id_o, access_time=access_data_s$access_time)


####���t�̏���####
##���t���f�[�g�^�ɕύX
regist_date <- as.Date(access_data_s$regist_date)
access_date <- as.Date(access_data_s$access_date)


##���Ԃ𐔒l�^�ɕύX���āA���ԑт���؂�
regist_time <- access_data_s$regist_time   #�o�^������unix�^�C��
access_time <- access_data_s$access_time   #�A�N�Z�X������unix�^�C��
access_hour_0 <- gsub(":", "", as.character(access_data_s$access_hour))
access_hour <- as.numeric(paste(gsub("0","", substring(access_hour_0, 1, 1)), substring(access_hour_0, 2), sep=""))

#���ԑт�5�ɋ�؂�
time_class <- ifelse(access_hour >= 0 & access_hour < 600, "����", ifelse(access_hour >= 600 & access_hour < 900, "��",
                                                                        ifelse(access_hour >= 900 & access_hour < 1800, 
                                                                               "�S��", "��")))
data.frame(time=access_data_s$access_hour, time_class, session_id)
table(time_class)   #���ԑт��ƂŃA�N�Z�X�����W�v


#���ԑт��_�~�[�ϐ���
time_table <- table(1:length(time_class), time_class)
time_name <- colnames(time_table)
time_d <- matrix(as.numeric(time_table), nrow=nrow(time_table), ncol=length(time_name))
colnames(time_d) <- time_name


##���t����j���_�~�[���쐬(�������x����)
weekday_o <- weekdays(access_date)
normalday <- c("���j��", "�Ηj��", "���j��", "�ؗj��")
holiday <- c("���j��", "�y�j��", "���j��")

#���������l�̃_�~�[
weekday <- ifelse(weekday_o %in% holiday, 1, 0)   


####�����ϐ��Ƒł��؂�w���ϐ��̍쐬####
user_unique <- unique(ID$user_id)   #���[�U�[id�̃��j�[�N
Y <- ifelse(access_data_s$vcs_view_limit_date == "", 0, 1)
T <- data.frame(id=as.character(user_id), Y, access_time)

#�ł��؂�x�N�g���̊i�[�p�p�����[�^
Z1 <- c()
Z2 <- c()
V <- c()

#���[�U�[���ƂɎw���ϐ����擾���āA�x�N�g���Ƃ��Č���
for(i in 1:length(user_unique)){
  print(i)
  T_ind <- T[T$id == user_unique[i], ]
  z1 <- rep(0, nrow(T_ind))
  z2 <- rep(0, nrow(T_ind))
  
  #�ł��؂�w���ϐ����쐬
  if(length(T_ind[T_ind$Y==1, 1])==0) {
    #y=1���Ȃ���΁A0�x�N�g���̂܂܌���
    Z1 <- c(Z1, z1)
    Z2 <- c(Z2, z2)
  } else {
    #y=1������΁A�ł��؂�ϐ��Ƃ��̌�A�N�Z�X�f�[�^�̎w���ϐ��𕪂���
    y_min <- min(T_ind[T_ind$Y==1, "access_time"])
    z1[subset(1:nrow(T_ind), T_ind$access_time==y_min)] <- 1
    if(which.max(z1)!=length(z1)) {z2[(which.max(z1)+1):length(z2)] <- 1}
    V <- rbind(V, c(i, length(z1), length(z2)))
    
    #�x�N�g��������
    Z1 <- c(Z1, z1)
    Z2 <- c(Z2, z2)
  }
}


##�h���C�����𑽒l�_�~�[�ϐ��ɕϊ�
domain <- access_data_s$from_domain
names(table(unique(domain)))
table(domain)


#�A�N�Z�X���̃C���f�b�N�X���쐬
index.facebook <- grep("facebook", domain)
index.google <- grep("google", domain)
index.yahoo <- grep("yahoo", domain)
index.bing <- grep("bing", domain)
index.docomo <- grep("docomo", domain)
index.au <- grep("auone", domain)
index.other <- c(index.facebook, index.google, index.yahoo, index.bing, index.docomo, index.au)


#�_�~�[�ϐ��ɕϊ�
#�p�����[�^�̊i�[�p�x�N�g��
facebook <- rep(0, length(domain))
google <- rep(0, length(domain))
yahoo <- rep(0, length(domain))
bing <- rep(0, length(domain))
docomo <- rep(0, length(domain))
au <- rep(0, length(domain))
other <- rep(0, length(domain))


#�h���C�����ʂ�01�_�~�[���i�[
facebook[index.facebook] <- 1
google[index.google] <- 1
yahoo[index.yahoo] <- 1
bing[index.bing] <- 1
docomo[index.docomo] <- 1
au[index.au] <- 1
other[-index.other] <- 1

#�f�[�^������
domain_d <- cbind(facebook, google, yahoo, bing, docomo, au, other)   
colSums(domain_d)

##���p�ړI���_�~�[�ϐ��ɕύX
unique(access_data_s$purpose)
purpose <- as.character(access_data_s$purpose)

#purpose�����l�������̂𕶎���ɖ߂�
pur_char <- ifelse(purpose=="10", "0010", ifelse(purpose=="100", "0100", ifelse(purpose=="1", "0001",
                                                                                ifelse(purpose=="101", "0101", ifelse(purpose=="11", "0011", ifelse(is.na(purpose)==TRUE, "0000",                                                                                                                                        purpose))))))


#������𕪊����ă_�~�[�ϐ�������
purpose_d <- matrix(0, nrow=length(purpose), ncol=4)
for(i in 1:4){
  purpose_d[, i] <- as.numeric(substr(pur_char, i, i))
}


####�y�[�W�����̃N�����W���O####
##�y�[�W�������_�~�[�ϐ���
no <- 1:nrow(access_data_s)
route <- access_data_s$route

#���l�_�~�[�ϐ���
route_cnt <- table(no, route)
route_name <- colnames(route_cnt)   #�y�[�W��
route_M <- matrix(as.numeric(route_cnt), nrow=nrow(route_cnt), ncol(route_cnt))   #�s��ɕϊ�
colnames(route_M) <- route_name


#kpi�֘A�y�[�W�̃_�~�[
answer <- c("answer", "survey_info")
teikei <- c("student_free", "tenshoku")
pay <- c("payment_info")

answer.v <- as.vector(rowSums(route_M[, answer]))
teikei.v <- as.vector(rowSums(route_M[, teikei]))
pay.v <- as.vector(route_M[, pay])
regist.v <- as.vector(route_M[, "regist_complete"])
pre.v <- as.vector(rowSums(route_M[, c(answer, teikei, pay)]))

#�A�N�Z�X����1000�����̃A�N�Z�X�̃y�[�W�͂��̑��ɕ���
kpi.gv <- setdiff(route_name, c(answer, teikei, pay))   #�s�v�ϐ�������
route_gv <- as.matrix(route_M[, kpi.gv])
gv_name <- colnames(route_gv)

#�A�N�Z�X����1000���ȏ�̃y�[�W�𒊏o
index.route <- subset(1:length(gv_name), colSums(route_gv) > 1000)   

#�A�N�Z�X����1000���ȏ�̃y�[�W�̃_�~�[
route1_v <- as.matrix(route_gv[, index.route])
route1_name <- colnames(route1_v)

#���̑��̃y�[�W�̃_�~�[
route2_v <- as.vector(rowSums(route_gv[, -index.route]))

#�f�[�^������
page_freq1 <- data.frame(answer=answer.v, teikei=teikei.v, pay=pay.v, regist=regist.v, route1_v, other=route2_v)
page_freq2 <- data.frame(pre.v, regist=regist.v, route1_v, other=route2_v)


####�{����Ƃ̃N�����W���O####
##��Ə�������
table(as.character(access_data_s$company_id1))
table(as.character(access_data_s$company_id2))

#��ƃA�N�Z�X�����f�[�^�t���[���ɕϊ�
company1 <- data.frame(no=1:nrow(access_data_s), master_id=access_data_s$company_id1)
company2 <- data.frame(no=1:nrow(access_data_s), act_id__c=access_data_s$company_id2)

#id�𕶎���ɕϊ�
company_vcs$act_id__c <- as.character(company_vcs$act_id__c)
company_vcs$master_id <- as.character(company_vcs$master_id)
company_all$master_id <- as.character(company_all$master_id)
company1$master_id <- as.character(company1$master_id)
company2$act_id__c <- as.character(company2$act_id__c)


#company_vcs��act_id��null�̍s���폜
index.null <- subset(1:nrow(company_vcs), company_vcs$act_id__c=="")   

##��Ə��Ɖ{����Ƃ�����
company1_join <- left_join(company1, company_all, by="master_id")
company2_join <- left_join(company2, company_vcs[-index.null, ], by="act_id__c")

#no�ŏd�������s���폜����
company1_join.dp <- company1_join[!duplicated(company1_join$no), ]
company2_join.dp <- company2_join[!duplicated(company2_join$no), ]

#���������f�[�^���m�F
company_j <- data.frame(id1=company1$master_id, id2=company2$act_id__c, route, mid1=company1_join.dp$master_id, 
                        mid2=company2_join.dp$master_id, aid=company2_join.dp$act_id__c, 
                        c1=company1_join.dp$company_name, c2=company2_join.dp$company_name)

#2��join����company��1�̍s��ɂ܂Ƃ߂�
index.join1 <- subset(1:nrow(company1_join.dp), company1_join.dp$master_id!="")
index.join2 <- subset(1:nrow(company2_join.dp), company2_join.dp$master_id!="")

#company�𓝍�
company_join <- company1_join.dp
company_join[index.join2, ] <- company2_join.dp[index.join2, -2]

company_master <- data.frame(id1=company1$master_id, id2=company2$act_id__c, route, company_join)

##���R�~�����ƂɃ_�~�[�ϐ���
vcs <- paste(company_join$public_vcs_flag, company_join$public_vcs_30_flag)
v_name <- c("���̑�", "30���ȏ�", "30������", "�Ȃ�")
vcs_d <- matrix(as.numeric(table(1:length(vcs), vcs)), nrow=length(vcs), ncol=length(v_name))
colnames(vcs_d) <- v_name
vcs_d <- vcs_d[, -1]
cbind(ID, vcs_d, Z1, Z2)


##�{���̎d�����擾
company_freq <- data.frame(id=user_id[Z2!=1], freq=company_join$master_id[Z2!=1])
company_table <- table(company_freq$id, company_freq$freq)
company_name <- colnames(company_table)
company_name <- company_name[-1]


#��Ɖ{���������s��
Cfreq_M0 <- matrix(as.numeric(company_table), nrow=nrow(company_table), ncol=ncol(company_table))
Cfreq_M <- Cfreq_M0[, -1]


#��Ɖ{����������{�����A�{����Ɛ��A�ő�{�������擾
company_sum <- rowSums(Cfreq_M)
company_max <- apply(Cfreq_M, 1, max)
company_cnt <- apply(Cfreq_M, 1, function(x) sum(x > 0))


#�f�[�^���m�F
#GLM�̃e�X�g
id_uni <- unique(ID$user_id)
Y_TEST <- c()

for(i in 1:length(id_uni)){
  yy <- max(Z1[ID$user_id==id_uni[i]])
  Y_TEST <- c(Y_TEST, yy)
  print(yy)
}

data.frame(ID$user_id, Z1)
DATA1 <- data.frame(Y=Y_TEST, sum=company_sum, max=company_max/company_sum, cnt=company_cnt/company_sum)
DATA2 <- data.frame(Y=Y_TEST, sum=company_sum, max=company_max, cnt=company_cnt)
TEST_DATA1 <- na.omit(DATA1)
TEST_DATA1 <- TEST_DATA1[TEST_DATA1$sum > 10, ]

res <- glm(Y ~ sum + cnt, data=TEST_DATA1, family=binomial)
summary(res)

##�Z�J���h�Z�b�V�����̎w���ϐ�
regist_time 
access_time


####�R���e���c�{�������̏���####
q_no <- as.character(access_data_s$q_no)
q_no[is.na(q_no)] <- 0   #na��0�ɂ��Ă���

#�_�~�[�ϐ���
q_no_table <- table(1:length(q_no), q_no)[, -1]
q_no_name <- colnames(q_no_table)
contents <- matrix(as.numeric(q_no_table), nrow=length(q_no), ncol=length(q_no_name))   #�s��ϊ�

#�f�[�^���m�F
colSums(contents)
q_no_name


##GLM�e�X�g2
X_fun <- as.matrix(data.frame(user_id, contents) %>% 
                     dplyr::group_by(user_id) %>%
                     dplyr::summarize_each(dplyr::funs(sum)))

DATA2 <- data.frame(Y=Y_TEST, sum=company_sum, max=company_max, cnt=company_cnt, c=X_fun[, -1])
res <- glm(Y ~ ., data=DATA2[, -3], family=binomial)
summary(res)


####�A���ϐ��̏���####
#�ݐσZ�b�V������
session_cnt <- session_id

#�o�^���̓o�^����
toroku <- log(access_data_s$toroku_time)
toroku[is.infinite(toroku)] <- 0


####��l�ϐ����_�~�[�ϐ��ɕϊ�####
keyword <- ifelse(access_data_s$keyword=="", 0, 1)   #�L�[���[�h�����������ǂ���  
sex <- ifelse(access_data_s$sex=="�j��", 1, 0)   #����
device <- ifelse(access_data_s$user_device=="sp", 1, 0)   #���[�U�[�f�o�C�X
career <- ifelse(access_data_s$career=="���K�Ј�", 1, 0)   #�ٗp�`��
mail_follow <- ifelse(access_data_s$mail_follow==1, 1, 0)   #���[���̃t�H���[�L��
j_change <- ifelse(is.na(access_data_s$job_change_timing), 0, 1)   #�]�E���������܂��Ă��邩�ǂ���


##���l�ϐ����_�~�[�ϐ��ɕϊ�
#�A�h�~�b�V�����^�C�v
admission <- matrix(as.numeric(table(1:length(user_id), access_data_s$admission_type)), nrow=length(user_id), ncol=3)


##�s���{�����O���[�v�����āA�_�~�[�ϐ���
pref <- access_data_s$pref
tokyo <- c("�����s")
kanto <- c("��ʌ�", "�_�ސ쌧", "��t��", "�Ȗ،�", "��錧", "�Q�n��")
kinki <- c("���{", "���s�{", "���Ɍ�", "�ޗǌ�", "���ꌧ", "�a�̎R��")
tyukyo <- c("���m��", "�O�d��", "�򕌌�")

#�s���{����n�悲�ƂɃO���[�v��
pref_o <- ifelse(pref %in% tokyo, "�����s",  ifelse(pref %in% kanto, "�֓���", ifelse(pref %in% kinki, "�ߋE��",
                                                                                ifelse(pref %in% tyukyo, "������", "���̑�"))))

#�s���{�����O���[�v�ʂɃ_�~�[�ϐ���
pref_table <- table(1:length(pref), pref_o)
pref_name <- colnames(pref_table)
pref_gr <- matrix(as.numeric(pref_table), nrow=length(pref), ncol=length(pref_name))   #�_�~�[�ϐ���
colnames(pref_gr) <- pref_name

colSums(pref_gr)   #�f�[�^���m�F


##�N���N��ɕϊ����ă_�~�[�ϐ���
age <- access_data_s$age 
age_class_o <- ifelse(age < 20, "10�Α�", ifelse(age >= 20 & age < 30, "20�Α�", ifelse(age >= 30 & age < 40, "30�Α�",
                                                                                    ifelse(age >= 40 & age < 50, "40�Α�", ifelse(age >= 50 & age < 60, "50�Α�", "60�Α�")))))
age_table <- table(1:length(age), age_class_o)
age_name <- colnames(age_table)

#�_�~�[�ϐ���
age_class <- matrix(as.numeric(age_table), nrow=length(age), ncol=length(age_name))
colnames(age_class) <- age_name


##�����ƊE���_�~�[�ϐ���
industry_o <- access_data_s$industry
industry_table <- table(1:length(industry_o), industry_o)
industry_name <- colnames(industry_table)

industry.m <- matrix(as.numeric(industry_table), nrow=length(industry_o), ncol=length(industry_name))   #�_�~�[�ϐ���
industry <- industry.m[, -1]
industry_name[-1]
colSums(industry)   #�f�[�^���m�F


##�E�Ƃ��_�~�[�ϐ���
occupation_o <- access_data_s$occupation_l
occupation_table <- table(1:length(occupation_o), occupation_o)
occupation_name <- colnames(occupation_table)

occupation.m <- matrix(as.numeric(occupation_table), nrow=length(occupation_o), ncol=length(occupation_name))   #�_�~�[�ϐ���
occupation <- occupation.m[, -1]
occupation_name[-1]
colSums(occupation)   #�f�[�^���m�F


##�N���N���X
income_o <- access_data_s$income_class
income_table <- table(1:length(income_o), income_o)
income_table <- income_table[, -1]
income_name <- colnames(income_table)

income <- matrix(as.numeric(income_table), nrow=length(income_o), ncol=length(income_name))   #�_�~�[�ϐ���
colSums(income)   #�f�[�^�̊m�F

##�f�o�C�X���p����
device_use <- matrix(0, nrow=n, ncol=3)


#���[�U�[���Ƃɗ��p�f�o�C�X���擾
for(i in 1:length(user_unique)){
  print(i)
  if(max(device[user_id==user_unique[i]])==1 & length(unique(device[user_id==user_unique[i]]))==1) {
    print("�X�}�z")
    device_use[i, 1] <- 1
    
  } else if(max(device[user_id==user_unique[i]])==0 & length(unique(device[user_id==user_unique[i]]))==1) {
    print("PC")
    device_use[i, 2] <- 1
    
  } else {
    print("����")
    device_use[i, 3] <- 1
  }
}
colnames(device_use) <- c("sp", "PC", "����")   #���O������
colSums(device_use)


####���v���f���̂��߂̃f�[�^���쐬####
####���ݕϐ����f�����O�̂��߂̐����ϐ����쐬####
##�f���O���t�B�b�N�ϐ�
Demo <- data.frame(sex, age_class, pref_gr, career, job=occupation, income=income, industry=industry, toroku=toroku,
                   ad=admission, mail=mail_follow, purpose=purpose_d, j_change=j_change)[Z2!=1, ]
demo_name <- colnames(Demo)


#���[�U�[���ƂɃf���O���t�B�b�N�ϐ���1�ɂ܂Ƃ߂�
id_z2 <- user_id[Z2!=1] 
Z.demo <- matrix(0, nrow=n, ncol=ncol(Demo))
len <- c()

for(i in 1:n){
  print(i)
  xz <- Demo[id_z2==user_unique[i], ]
  len <- c(len, nrow(xz))
  Z.demo[i, ] <- as.numeric(xz[1, ])
}

colnames(Z.demo) <- demo_name


#�����ϐ���ݒ�
Y.demo <- data.frame(Y, ID) %>%
  dplyr::group_by(user_id) %>%
  dplyr::summarize(Y=max(Y))
Yi <- as.matrix(Y.demo)[, 2]

####�������f���̂��߂̃p�l���f�[�^���쐬####
##���[�U�[�o�^�O�Ƃ��Ƃ��s�ʂ���
regist_flag <- c()
for(i in 1:n){
  print(i)
  
  #user_id���Ƃɓo�^�O��̃t���b�O�����Ă�
  index.ind <- subset(1:length(route), ID[, 2] == user_unique[i])
  y_regist <- rep(0, length(index.ind))
  r_comp <- which.max(route[index.ind]=="regist_complete")
  route[index.ind]
  
  if(r_comp > 1) {
    y_regist[r_comp:length(index.ind)] <- 1
  } else {
    diff <- access_time[index.ind] - regist_time[index.ind]
    y_regist[1:length(y_regist)] <- ifelse(diff >= -0.01, 1, 0)
  }
  #�x�N�g��������
  regist_flag <- c(regist_flag, y_regist)
}

##���[�U�[�o�^����24���Ԉȓ��Ƃ���ȏ�ɕ�����
diff <- access_time - regist_time   #�o�^���ԂƃA�N�Z�X���Ԃ̍���
day_s <- 60*60*24   #1����b���Z
SS_flag <- ifelse(diff >= day_s, 1, 0)


##���[�U�[�o�^�O�̍s�����L�^
#�t�@�[�X�g�����f�B�����烆�[�U�[�o�^�܂ł̌o�ߓ���
#���[�U�[�o�^�O��cv�֘A�y�[�W�𓥂�ł��邩�ǂ���
keika_time <- c()
keika_ind <- c()
kpi_cnt <- c()
kpi_pv <- ifelse(route %in% c(answer, teikei, pay), 1, 0)   #�S�̃A�N�Z�X��kpi�A�N�Z�X�L��

for(i in 1:n){
  print(i)
  
  #���[�U�[id�ʂɃt�@�[�X�g�����f�B���O����o�^�܂ł̌o�ߎ��Ԃ��L�^
  index.ind <- subset(1:length(access_time), ID[, 2] == user_unique[i])
  keika <- rep(0, length(index.ind))
  p_time <- abs(min(access_time[user_id==user_unique[i]] - regist_time[user_id==user_unique[i]])/60/60/24)
  keika[1:length(keika)] <- p_time
  keika_ind <- c(keika_ind, p_time)
  keika_time <- c(keika_time, keika)   #�x�N�g���Ɍ���
  
  #���[�U�[�o�^�O��cv�֘A�y�[�W�𓥂�ł��邩�ǂ������L�^
  k_cnt <- kpi_pv[user_id==user_unique[i] & regist_flag==0]
  kpi_cnt <- c(kpi_cnt, sum(k_cnt))   #�x�N�g���Ɍ���
}

##���[�U�[�o�^�O�̑�pv�����L�^
unregist_pv <- table(user_id, regist_flag)
ur_pv <- matrix(as.numeric(table(user_id, regist_flag)), nrow=length(user_unique), ncol=2)[, 1]


##���[�U�[�o�^�O�̃R���e���c�̉{����
q_cnt <- table(user_id[regist_flag==0], q_no[regist_flag==0])
index_q <- subset(1:n, 1:n %in% as.numeric(rownames(q_cnt)) == FALSE)

#�������Ă���id��0�x�N�g������
cont_q <- matrix(0, nrow=n, ncol=ncol(q_cnt))

cont_q[-index_q, ] <- q_cnt
cont_q[index_q, ] <- rep(0, ncol(q_cnt))
cont_m <- cont_q[, -1]

##���[�U�[�o�^�O�̃y�[�W�{����
pagename <- colnames(page_freq1)
page_unit <- as.numeric(as.matrix(page_freq1) %*% 1:length(pagename))

#���[�U�[�o�^�O�̃y�[�W�̉{����
p_cnt <- table(c(user_id[regist_flag==0], rep(3001, ncol(page_freq1))), c(page_unit[regist_flag==0], 1:length(pagename)))
p_cnt <- p_cnt[-nrow(p_cnt), ]
index_p <- subset(1:n, 1:n %in% as.numeric(rownames(p_cnt)) == FALSE)

#�������Ă���id��0�x�N�g������
page_q <- matrix(0, nrow=n, ncol(page_freq1))

page_q[-index_p, ] <- p_cnt
page_q[index_p, ] <- rep(0, ncol(page_freq1))
colnames(page_q) <- pagename


##�Z�b�V����id���ƂɃf�[�^���܂Ƃ߂ăp�l���f�[�^�ɂ���
#���[�U�[ID�A�Z�b�V�������ƂɃC���f�b�N�X���L�^����
session.list <- list()

#id�A�Z�b�V�������ƂɃ��X�g��
for(i in 1:n){
  print(i)
  unique_session <- unique(session_id[user_id==i])
  session.ind <- list()
  
  for(j in 1:length(unique_session)){
    index.session <- subset(1:length(session_id), user_id==i & session_id==j)
    session.ind[[j]] <- index.session
  }
  session.list[[i]] <- session.ind
}


#�Z�b�V����ID��1�s�ɂ܂Ƃ߂�ID���쐬
ID_D <- as.matrix(ID %>% 
                    dplyr::group_by(user_id) %>%
                    dplyr::count(session_id))[, -3]

ID.v <- data.frame(no=1:nrow(ID_D), ID_D)


##�Z�b�V����ID���ƂɃA�N�Z�X���A�A�N�Z�X�T�A�A�N�Z�X���ԂȂǃZ�b�V�����ň�ӂȕϐ����擾
#���o����Ώۂ̃f�[�^
TIME_DATA <- data.frame(ID, Y, Z1, Z2, regi=regist_flag, access_date, weekday, time_d, SS=SS_flag, 
                        domain=domain_d, q=contents, p=page_freq1, s_cnt=session_cnt-1, device)

#ID�A�Z�b�V�������Ƃɐ����ϐ����L�^
XZ <- list()
Cont <- list()
Page <- list()
LEN <- 0

for(i in 1:n){
  print(i)
  len <- length(session.list[[i]])
  for(j in 1:len){
    LEN <- LEN + 1
    Page[[LEN]] <- colSums(rbind(page_freq1[unlist(session.list[[i]][j]), ], rep(0, ncol(page_freq1))))
    Cont[[LEN]] <- colSums(rbind(contents[unlist(session.list[[i]][j]), ], rep(0, ncol(contents))))
    
    ind.session <- TIME_DATA[unlist(session.list[[i]][j]), ]
    index.ind <- subset(1:length(ind.session$Z1), ind.session$Z1==1)
    if(length(index.ind) > 0) {
      XZ[[LEN]] <- TIME_DATA[unlist(session.list[[i]][j])[index.ind[1]], ]
    } else {
      XZ[[LEN]] <- TIME_DATA[unlist(session.list[[i]][j])[1], ]
    }
  }
}

#���X�g�`�����s��`���ɕύX
X.panel1 <- do.call(rbind, XZ)
C.panel1 <- do.call(rbind, Cont)
P.panel1 <- do.call(rbind, Page)


##�ړI�A��Ɖ{�������̃f�[�^��R�t����
#��Ɖ{��������R�t����
com <- matrix(0, nrow=nrow(X.panel1), ncol=3)
company_cm <- cbind(company_cnt, company_max, company_sum)

for(i in 1:n){
  index.co <- subset(1:nrow(com), ID.v$user_id==i)
  com[index.co, ] <- matrix(company_cm[i, ], nrow=length(index.co), ncol=3, byrow=TRUE)
}


#�ړI��R�t����
pur <- matrix(0, nrow=nrow(X.panel1), ncol=ncol(purpose_d))

for(i in 1:n){
  index.pur <- subset(1:nrow(pur), ID.v$user_id==i)
  if(length(index.pur) > 1) {
    pur[index.pur, ] <-  matrix(purpose_d[ID$user_id==i, ][1, ], nrow=length(index.pur), ncol=ncol(purpose_d), byrow=TRUE)
  } else {
    pur[index.pur, ] <-  matrix(purpose_d[ID$user_id==i, ], nrow=length(index.pur), ncol=ncol(purpose_d), byrow=TRUE)
  } 
}


##�y�[�W�{�������A�R���e���c�{���������L�^����
C.id1 <- cbind(X.panel1[, c("user_id", "session_id")], cn=C.panel1, cp=matrix(0, nrow(C.panel1), ncol(C.panel1)))
C.id2 <- cbind(C.id1, cf=matrix(0, nrow(C.panel1), ncol(C.panel1)))
P.id1 <- cbind(X.panel1[, c("user_id", "session_id")], n=P.panel1, p=matrix(0, nrow(P.panel1), ncol(P.panel1)))
P.id2 <- cbind(P.id1, pf=matrix(0, nrow(P.panel1), ncol(P.panel1)))

r1 <- 3:((ncol(P.id1)-2)/2+2)
r2 <- ((ncol(P.id1)-2)/2+3):ncol(P.id1)
r3 <- (ncol(P.id1)+1):ncol(P.id2)

for(i in 1:n){  
  print(i)
  #�ߋ��̉{����������
  cont_hist <- cont_m[i, ]
  page_hist <- page_q[i, ]
  index.user <- subset(1:nrow(C.id1), C.id1$user_id==i)
  user_session <- length(index.user)
  
  for(j in 1:user_session){
    C.id1[index.user[j], 13:ncol(C.id1)] <- cont_hist
    C.id2[index.user[j], 13:ncol(C.id1)] <- cont_hist
    P.id1[index.user[j], r2] <- page_hist
    P.id2[index.user[j], r2] <- page_hist
    
    #�ߋ��̉{�������ƌ��݂̉{���������擾
    c.past <- C.id1[index.user[j], 13:ncol(C.id1)] 
    c.now <- C.id1[index.user[j], 3:12]
    p.past <- P.id1[index.user[j], r2] 
    p.now <- P.id1[index.user[j], r1]
    
    #���߂Č���Ȃ�p�x�𐔂���
    c.first_cnt <- abs(c.past - c.now)
    c.first_flag <- c.past-c.now < 0 & c.past==0
    c.first_cnt[c.first_flag==FALSE]  <- 0 
    
    p.first_cnt <- abs(p.past - p.now)
    p.first_flag <- p.past-p.now < 0 & p.past==0
    p.first_cnt[p.first_flag==FALSE]  <- 0
    
    #����{���f�[�^����
    C.id2[index.user[j], (ncol(C.id1)+1):ncol(C.id2)] <- c.first_cnt
    P.id2[index.user[j], r3] <- p.first_cnt
    
    #�{���������X�V
    if(j==user_session) {
      break
    } else {
      cont_hist <- C.id2[index.user[j], 13:ncol(C.id1)] + C.id2[index.user[j], 3:12]
      page_hist <- P.id2[index.user[j], r2] + P.id2[index.user[j], r1]
    } 
  }
}


##�ϐ��̃X�P�[����
cscale <- scale(C.id2[, 13:22]) 
cp <- cscale + matrix(abs(apply(cscale, 2, min)), nrow=nrow(cscale), ncol=ncol(cscale), byrow=TRUE)
c.ones <- ifelse(cp > 0, 1, 0)   #�ߋ��ɉ{�����Ă��邩�ǂ���

pscale <- scale(P.id2[, r2])
pp <- pscale + matrix(abs(apply(pscale, 2, min)), nrow=nrow(pscale), ncol=ncol(pscale), byrow=TRUE)
p.ones <- ifelse(P.id2[, r2] > 0, 1, 0)

p_cnt.scale <- scale(rowSums(P.id2[, 33:62]))
p_cnt <- p_cnt.scale + abs(min(p_cnt.scale))

com.rate <- com[, 1:2]/com[, 3]
com.rate[is.nan(com.rate)] <- 0
com_cnt <- scale(com[, 3]) + abs(min(scale(com[, 3])))


####�������f���ƃN���X�Z�N�V�������f���𓖂Ă͂߂�####
##�p�l���f�[�^���쐬���Đ������f���𓖂Ă͂߂�
#�f�[�^�̐ݒ�
ID.v1 <- X.panel1[, 2:3]   #ID�̐ݒ�

fit_data1 <- data.frame(Y=X.panel1$Z1, X.panel1[, c(3, 9:21)], C.id2[, 3:12], C.id2[, 23:32],   #�Ó� 
                        P.id2[, c(69, 80, 81, 87)], P.id2[, c(9, 20, 21, 27)], pur=pur, cnt=p_cnt)

fit_data2 <- data.frame(Y=X.panel1$Z1, X.panel1[, c(3, 9:21)], C.id2[, 23:32], ones=c.ones,
                        P.id2[, c(69, 80, 81, 87)], P.id2[, c(9, 20, 21, 27)], pur=pur, cnt=p_cnt)

fit_data3 <- data.frame(Y=X.panel1$Z1, X.panel1[, c(3, 9:21)], C.id2[, 3:12], C.id2[, 23:32], ones=c.ones,
                        P.id2[, c(69, 80, 81, 87)], P.id2[, c(9, 20, 21, 27)], pur=pur, cnt=p_cnt)


ID.v2 <- ID.v1[X.panel1$Z2==0, ]
fit_data01 <- fit_data1[X.panel1$Z2==0, -c(5, 15)]
fit_data02 <- fit_data2[X.panel1$Z2==0, -c(5, 15)]
fit_data03 <- fit_data3[X.panel1$Z2==0, -c(5, 15)]

#GLM�𓖂Ă͂߂�
fit1 <- glm(Y ~ ., data=fit_data01, family=binomial)
fit2 <- glm(Y ~ ., data=fit_data02, family=binomial)
fit3 <- glm(Y ~ ., data=fit_data03, family=binomial)

summary(fit1)
summary(fit2)
summary(fit3)


##�f���O���t�B�b�N�ϐ���GLM�𓖂Ă͂߂�(�\������)
#�f�[�^�̐ݒ�
Z.demo_add <- data.frame(Z.demo, cnt=company_cnt, sum=company_sum, device_use)

#�璷�ȕϐ�������
index.z <- setdiff(colnames(Z.demo_add), c("X10�Α�", "���̑�", "job.12", "ad.2", "income.6", "industry.9", "PC"))   
Zi.demo <- Z.demo_add[, index.z]

res <- glm(Y ~ ., data = data.frame(Y=Yi, Zi.demo[, -c(12:25)]), family=binomial)
summary(res)


####EM�A���S���Y����Split_Hazard_model�𐄒�####
####EM�A���S���Y���ŗ��p����֐��̐ݒ�####
##���W�b�g���f���̑ΐ��ޓx�̒�`
logit_LL <- function(x, Y, X){
  #�p�����[�^�̐ݒ�
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #���p�֐��̒�`
  U <- b0 + as.matrix(X) %*% b1
  
  #�ΐ��ޓx���v�Z
  Pr <- exp(U)/(1 + exp(U))   #�m���̌v�Z
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

##���S�f�[�^�ł�Split_hazard_mode�̑ΐ��ޓx
#�p�����[�^�̐ݒ�
cll <- function(b, Y, Yn, X, zpt, zk){
  #�p�����[�^�̐ݒ�
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #���p�֐����`
  U <- beta0 + as.matrix(X) %*% beta1
  
  #���W�X�e�B�b�N���f���̊m���Ƒΐ��ޓx�̌v�Z
  Pr <- exp(U)/(1 + exp(U))   #�m���̌v�Z
  LLc <- Y*log(Pr) + (1-Y)*log(1-Pr)
  
  #�[���ޓx�̌v�Z
  LLzero <- dbinom((1-Yn), 1, 1)   #�[���ޓx�̌v�Z
  
  LL <- sum(zpt[, 1]*LLc)
  return(LL)
}


##�ϑ��f�[�^�ł̖ޓx�Ɛ��ݕϐ�z�̌v�Z
ollz <- function(b, Y, Yn, X, r, id, n, zk){
  
  #�p�����[�^�̐ݒ�
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #�Z�O�����g���Ƃ̌��p�֐����`
  U <- beta0 + as.matrix(X) %*% beta1
  
  #���W�X�e�B�b�N���f���̊m���Ƒΐ��ޓx�̌v�Z
  Pr <- exp(U)/(1 + exp(U))   #�m���̌v�Z
  LLi <- Pr^Y * (1-Pr)^(1-Y)   #�ޓx�̌v�Z
  LLzero <- dbinom((1-Yn), 1, 1)   #�[���ޓx�̌v�Z
  LCo <- cbind(LLi, LLzero)   #�ޓx�̌���
  
  #ID�ʂɖޓx�̐ς����
  LLho <- matrix(0, nrow=n, ncol=zk)
  for(i in 1:n){
    if(sum(id==i)==1){
      LLho[i, ] <- LCo[id==i, ]
    } else {
      LLho[i, ] <- apply(LCo[id==i, ], 2, prod) 
    }
  }
  
  #�ϑ��f�[�^�ł̑ΐ��ޓx
  LLo <- sum(log(apply(r * LLho, 1, sum)))
  
  #���ݕϐ�z�̌v�Z
  z0 <- r * LLho   #���ݕϐ�z�̕��q
  z1 <- z0 / rowSums(z0)   #���ݕϐ�z�̌v�Z
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}

####EM�A���S���Y���̐ݒ�Ə���####
##EM�A���S���Y���̏����l�̐ݒ�
#�f�[�^�̐ݒ�
ID.data <- data.frame(no=1:nrow(ID.v2), ID.v2)   #ID�ƃZ�b�V�����ԍ�
id <- ID.data$user_id
Y1 <- fit_data01[, 1]   #�������f���̉����ϐ�
XM <- fit_data01[, -c(1:3)]   #�������f���̐����ϐ�
XM[, 11:20] <- XM[, 11:20] - XM[, 21:30]
XM[, 35:38] <- XM[, 35:38] - XM[, 31:34]

Zi <- Zi.demo   #���݃��f���̐����ϐ�
Y2 <- Yi
zpt <- matrix(0, nrow=nrow(XM), ncol=2)

#��������
#Zi1 <- Zi
Zi[is.na(Zi)] <- 0
XM[is.na(XM)] <- 0
#Zi <- Zi[, -(34:42)]

#�[���ޓx�̂��߂�y�̐ݒ�
index.id <- subset(ID.data$user_id, Y1==1)
Y.zeros <- ifelse(ID.data$user_id %in% index.id, 1, 0)


####EM�A���S���Y���ŃX�v���b�g�n�U�[�h���f���𐄒�####
##�ŗǒl���o��܂Ŕ������J��Ԃ�
rt <- 40   #�J��Ԃ���
LL.process <- c()
LL.list <- c()
beta.list <- matrix(0, nrow=rt, ncol=ncol(XM)+1)
theta.list <- matrix(0, nrow=rt, ncol=ncol(Zi)+1)
z.list <- array(0, dim=c(n, 2, rt))

for(rp in 1:rt){
  print(rp)
  
  ##EM�A���S���Y���̏����l��ݒ�
  #�������f���̏����l��ݒ�
  res.surv <- glm(Y1 ~ ., data=XM, family=binomial(link="logit"))
  summary(res.surv)
  beta <- coef(res.surv)
  
  #���݃��f���̏����l��ݒ�
  res.latent <- glm(Y2 ~ ., data=Zi, family=binomial(link="logit"))
  summary(res.latent)
  theta <- coef(res.latent)
  
  #�m���̌v�Z
  Z <- as.matrix(cbind(1, Zi)) %*% theta   #���W�b�g�̌v�Z
  Pr.z <- exp(Z)/(1+exp(Z))   #�m���̌v�Z
  r <- cbind(Pr.z, 1-Pr.z)   #�������̌v�Z
  
  ##1��ڂ̃X�e�b�v�̂�M�X�e�b�v�������l�Ɉˑ����Ȃ�SANN�Ő���
  #���ݕϐ�z�Ɗϑ��f�[�^�̑ΐ��ޓx�̏����l�̐ݒ�
  
  oll <- ollz(b=beta, Y=Y1, Yn=Y.zeros, X=as.matrix(XM), r=r, id=ID.data$user_id, n=n, zk=2)
  z <- oll$z1
  LL1 <- oll$LLo
  
  #���ݕϐ�z���p�l���`���ɕύX
  for(i in 1:n){
    index.id <- subset(1:length(id), id==i)
    zpt[index.id, ] <- matrix(z[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
  }
  
  
  #���S�f�[�^�ł̐������f���̍Ŗސ���(M�X�e�b�v)
  beta.first <- beta   #�����l��ݒ�
  
  for(i in 1:10000){
    beta1 <- beta.first + runif(length(beta), -0.7, 0.7)
    res <- try(optim(beta1, cll, Y=Y1, Yn=Y.zeros, X=XM, zpt=zpt, zk=2, method="Nelder-Mead", 
                     hessian=FALSE, control=list(fnscale=-1, trace=TRUE, maxit=3000)), silent=TRUE)
    if(class(res) == "try-error") {next} else {break}   #�G���[����
  }
  beta <- res$par
  
  #������r�̍X�V
  res.latent <- glm(z ~ ., data=Zi, family=binomial(link="logit"))
  theta <- coef(res.latent)    #���ݕϐ��̃p�����[�^
  
  #�m���̌v�Z
  Z <- as.matrix(cbind(1, Zi)) %*% theta   #���W�b�g�̌v�Z
  Pr.z <- exp(Z)/(1+exp(Z))   #�m���̌v�Z
  r <- cbind(Pr.z, 1-Pr.z)   #�������̌v�Z
  
  #���ݕϐ�z�Ɗϑ��f�[�^�̑ΐ��ޓx�̏����l�̐ݒ�
  oll <- ollz(b=beta, Y=Y1, Yn=Y.zeros, X=as.matrix(XM), r=r, id=ID.data$user_id, n=n, zk=2)
  z <- oll$z1
  LL1 <- oll$LLo
  
  
  #EM�A���S���Y���̐ݒ�
  iter <- 0
  dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l��ݒ�
  tol <- 1   
  
  ##EM�A���S���Y���ŃX�v���b�g�n�U�[�h���f���𐄒�
  while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
    #���ݕϐ�z���p�l���`���ɕύX
    for(i in 1:n){
      index.id <- subset(1:length(id), id==i)
      zpt[index.id, ] <- matrix(z[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
    }
    
    #���S�f�[�^�ł̐������f���̍Ŗސ���(M�X�e�b�v)
    res <- optim(beta, cll, Y=Y1, Yn=Y.zeros, X=XM, zpt=zpt, zk=2, method="Nelder-Mead", 
                 hessian=FALSE, control=list(fnscale=-1, maxit=2000))
    beta <- res$par
    
    #������r�̍X�V
    res.latent <- glm(z ~ ., data=Zi, family=binomial(link="logit"))
    theta <- coef(res.latent)    #���ݕϐ��̃p�����[�^
    Z <- as.matrix(cbind(1, Zi)) %*% theta   #���W�b�g�̌v�Z
    Pr.z <- exp(Z)/(1+exp(Z))   #�m���̌v�Z
    r <- cbind(Pr.z, 1-Pr.z)   #�������̌v�Z
    
    #���ݕϐ�z�Ɗϑ��f�[�^�̑ΐ��ޓx���X�V
    oll <- ollz(b=beta, Y=Y1, Yn=Y.zeros, X=as.matrix(XM), r=r, id=ID.data$user_id, n=n, zk=2)
    z <- oll$z1
    LL <- oll$LLo
    
    #EM�A���S���Y���̃p�����[�^�̍X�V
    iter <- iter+1
    dl <- LL-LL1
    LL1 <- LL
    LL.process <- c(LL.process, LL)
    print(LL)
    if(abs(dl) < 5 & LL < max(LL.process)) {break}
  }
  
  #�p�����[�^����
  LL.list <- c(LL.list, LL) 
  beta.list[rp, ] <- beta
  theta.list[rp, ] <- theta
  z.list[, , rp] <- z 
}

beta <- res$par

##�œK�ȃp�����[�^������
LL.list
opt.LL <- which.max(LL.list)
LL.list[opt.LL]   #��ă��f���̍ő�ΐ��ޓx
logLik(res.surv)   #���W�X�e�B�b�N��A�̑ΐ��ޓx

#���肳�ꂽ�p�����[�^
beta.best <- beta.list[opt.LL, ]
theta.best <- theta.list[opt.LL, ]
z.best <- z.list[, , opt.LL]


#��A�W���ɖ��O������
beta.name <- c("�ؕ�", colnames(XM))
theta.name <- c("�ؕ�", colnames(Zi))
names(beta.best) <- beta.name
names(theta.best) <- theta.name

#��A�W�����m�F
round(beta.best, 3)
round(coef(res.surv), 3)
round(theta.best, 3)

occupation_name

colSums(Zi)

round(zzz <- data.frame(pre=z.best[, 1], non_pre=z.best[, 2], Y2), 3)   #���ݕϐ��Ɣ����ϐ��̔�r
colSums(z.best)/n   #������

colSums(ifelse(X[, 11:30] > 0, 1, 0))
summary(XM)

##�K���x��t����
round(res.surv$aic, 2)
round(AIC <- -2*LL.list[opt.LL] + 2*length(beta.best), 2)   ##AIC


res.latent <- glm(z.best ~ ., data=Zi, family=binomial(link="logit"))
summary(res.latent)


for(i in 1:n){
  index.id <- subset(1:length(id), id==i)
  zpt[index.id, ] <- matrix(z.best[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
}

#���S�f�[�^�ł̐������f���̍Ŗސ���(M�X�e�b�v)
res <- optim(beta.best, cll, Y=Y1, Yn=Y.zeros, X=XM, zpt=zpt, zk=2, method="CG", 
             hessian=TRUE, control=list(fnscale=-1, maxit=2000))
res$par
beta.best


apply(Zi, 2, quantile)

colSums(Zi)
apply(XM[, ], 2, quantile)
colSums(XM)

test_data <- XM[, 31:36]

uni.id <- unique(id)
sums <- matrix(0, nrow=length(uni.id), ncol=ncol(XM))

for(i in 1:length(uni.id)){
  print(i)
  sums[i, ] <- colSums(XM[id==i, ] > 0)
}
colnames(sums) <- colnames(XM)

colSums(Zi)

nrow(XM)-1893-8388-3719-164-82-33

as.matrix(XM) %*% beta.best


for(i in 1:n){
  index.id <- subset(1:length(id), id==i)
  zpt[index.id, ] <- matrix(z.best[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
}


apply(XM, 2, quantile)
colSums(XM > 0)

Zii <- data.frame(id=1:nrow(Zi), Zi)

XZ <- left_join(data.frame(id, XM), Zii, by="id")
XZ[0, ]

logit <- as.matrix(cbind(1, XM)) %*% beta.best
XZM <- round(cbind(no=1:length(id), id, session=X.panel1[X.panel1$Z2==0, "session_id"], Y1, p1=exp(logit)/(1+exp(logit)), 
                   p2=res.surv$fitted.values, z=zpt, XZ[, -1]), 2)

path_output <- paste("C:/Users/okuno/Desktop/", "�J���e", ".csv", sep="")   #�ۑ��p�X
write.table(XZM, path_output, row.names=FALSE, col.names = TRUE,
            append=FALSE, sep=",", fileEncoding="Shift-JIS")   #�����o��


length(res.surv$fitted.values)
length(Y1)
nrow(XM)



