#####カルテ利用性向上のための潜在変数モデリング#####
library(MASS)
library(survival)
library(dplyr)
library(caret)
library(ggplot2)
library(lattice)

set.seed(3181)

####データの読み込みとデータをサンプリング####
access_all <- read.csv("C:/Statistics/dat/access_log.csv", fileEncoding = "Shift-JIS")
company_all <- read.csv("C:/Statistics/dat/company_without_act_id.csv", fileEncoding = "Shift-JIS")
company_vcs <- read.csv("C:/Statistics/dat/company_with_act_id.csv", fileEncoding = "Shift-JIS")
access_all1 <- access_all

index.user1 <- subset(names(sort(table(access_all1$user_id), decreasing=TRUE)), 
                      sort(table(access_all1$user_id), decreasing=TRUE) < 500)
access_all <- access_all1[access_all1$user_id %in% index.user1, ]
sort(table(access_all$user_id), decreasing=TRUE)


#不要変数を除去
gv <- setdiff(colnames(access_all), c("u_regist_date", "url", "user_agent", "purpose1", "purpose2",
                                      "purpose3", "purpose4", "purpose5", "purpose6"))   #不要変数を除去


access_data <- access_all[, gv]


##ユーザーidをランダムサンプリング(5000人)
n <- 10000
id_list <- unique(access_data$user_id)
id_sample <- sample(id_list, n)

#ランダムサンプリングを実行
access_data_s <- access_data[access_data$user_id %in% id_sample, ]


####データのクレンジング####
####IDの作成####
##user_idを数値に変更
td_id <- as.numeric(access_data_s$td_client_id)
user_id_o <- as.character(access_data_s$user_id)
unique_id <- unique(user_id_o)

#user_idを連番に変更
user_id <- array(0, dim=c(length(user_id_o)))
for(i in 1:length(unique_id)){
  print(i)
  user_id[user_id_o==unique_id[i]] <- i
}

##session_idを数値に変更
session_id_o <- as.numeric(access_data_s$session_id)
unique_id <- unique(user_id)

#session_idを連番に変更
session_id <- array(0, dim=c(length(session_id_o)))
for(i in 1:length(unique_id)){
  print(i)
  unique_session <- unique(session_id_o[user_id==unique_id[i]])
  
  #ユーザーidごとにセッションidの連番をつける
  for(j in 1:length(unique_session)){
    session_id[session_id_o==unique_session[j]] <- j
  }
}


#IDの結合と結合したIDの確認
ID <- data.frame(no=1:length(user_id), user_id, session_id)
ID_watch <- data.frame(user_id, session_id, user_id_o, session_id_o, access_time=access_data_s$access_time)


####日付の処理####
##日付をデート型に変更
regist_date <- as.Date(access_data_s$regist_date)
access_date <- as.Date(access_data_s$access_date)


##時間を数値型に変更して、時間帯を区切る
regist_time <- access_data_s$regist_time   #登録時期のunixタイム
access_time <- access_data_s$access_time   #アクセス時期のunixタイム
access_hour_0 <- gsub(":", "", as.character(access_data_s$access_hour))
access_hour <- as.numeric(paste(gsub("0","", substring(access_hour_0, 1, 1)), substring(access_hour_0, 2), sep=""))

#時間帯を5つに区切る
time_class <- ifelse(access_hour >= 0 & access_hour < 600, "未明", ifelse(access_hour >= 600 & access_hour < 900, "朝",
                                                                        ifelse(access_hour >= 900 & access_hour < 1800, 
                                                                               "全日", "夜")))
data.frame(time=access_data_s$access_hour, time_class, session_id)
table(time_class)   #時間帯ごとでアクセス数を集計


#時間帯をダミー変数化
time_table <- table(1:length(time_class), time_class)
time_name <- colnames(time_table)
time_d <- matrix(as.numeric(time_table), nrow=nrow(time_table), ncol=length(time_name))
colnames(time_d) <- time_name


##日付から曜日ダミーを作成(平日か休日か)
weekday_o <- weekdays(access_date)
normalday <- c("月曜日", "火曜日", "水曜日", "木曜日")
holiday <- c("金曜日", "土曜日", "日曜日")

#平日か求人のダミー
weekday <- ifelse(weekday_o %in% holiday, 1, 0)   


####応答変数と打ち切り指示変数の作成####
user_unique <- unique(ID$user_id)   #ユーザーidのユニーク
Y <- ifelse(access_data_s$vcs_view_limit_date == "", 0, 1)
T <- data.frame(id=as.character(user_id), Y, access_time)

#打ち切りベクトルの格納用パラメータ
Z1 <- c()
Z2 <- c()
V <- c()

#ユーザーごとに指示変数を取得して、ベクトルとして結合
for(i in 1:length(user_unique)){
  print(i)
  T_ind <- T[T$id == user_unique[i], ]
  z1 <- rep(0, nrow(T_ind))
  z2 <- rep(0, nrow(T_ind))
  
  #打ち切り指示変数を作成
  if(length(T_ind[T_ind$Y==1, 1])==0) {
    #y=1がなければ、0ベクトルのまま結合
    Z1 <- c(Z1, z1)
    Z2 <- c(Z2, z2)
  } else {
    #y=1があれば、打ち切り変数とその後アクセスデータの指示変数を分ける
    y_min <- min(T_ind[T_ind$Y==1, "access_time"])
    z1[subset(1:nrow(T_ind), T_ind$access_time==y_min)] <- 1
    if(which.max(z1)!=length(z1)) {z2[(which.max(z1)+1):length(z2)] <- 1}
    V <- rbind(V, c(i, length(z1), length(z2)))
    
    #ベクトルを結合
    Z1 <- c(Z1, z1)
    Z2 <- c(Z2, z2)
  }
}


##ドメイン元を多値ダミー変数に変換
domain <- access_data_s$from_domain
names(table(unique(domain)))
table(domain)


#アクセス元のインデックスを作成
index.facebook <- grep("facebook", domain)
index.google <- grep("google", domain)
index.yahoo <- grep("yahoo", domain)
index.bing <- grep("bing", domain)
index.docomo <- grep("docomo", domain)
index.au <- grep("auone", domain)
index.other <- c(index.facebook, index.google, index.yahoo, index.bing, index.docomo, index.au)


#ダミー変数に変換
#パラメータの格納用ベクトル
facebook <- rep(0, length(domain))
google <- rep(0, length(domain))
yahoo <- rep(0, length(domain))
bing <- rep(0, length(domain))
docomo <- rep(0, length(domain))
au <- rep(0, length(domain))
other <- rep(0, length(domain))


#ドメイン元別に01ダミーを格納
facebook[index.facebook] <- 1
google[index.google] <- 1
yahoo[index.yahoo] <- 1
bing[index.bing] <- 1
docomo[index.docomo] <- 1
au[index.au] <- 1
other[-index.other] <- 1

#データを結合
domain_d <- cbind(facebook, google, yahoo, bing, docomo, au, other)   
colSums(domain_d)

##利用目的をダミー変数に変更
unique(access_data_s$purpose)
purpose <- as.character(access_data_s$purpose)

#purposeが数値化したのを文字列に戻す
pur_char <- ifelse(purpose=="10", "0010", ifelse(purpose=="100", "0100", ifelse(purpose=="1", "0001",
                                                                                ifelse(purpose=="101", "0101", ifelse(purpose=="11", "0011", ifelse(is.na(purpose)==TRUE, "0000",                                                                                                                                        purpose))))))


#文字列を分割してダミー変数化する
purpose_d <- matrix(0, nrow=length(purpose), ncol=4)
for(i in 1:4){
  purpose_d[, i] <- as.numeric(substr(pur_char, i, i))
}


####ページ履歴のクレンジング####
##ページ履歴をダミー変数化
no <- 1:nrow(access_data_s)
route <- access_data_s$route

#多値ダミー変数化
route_cnt <- table(no, route)
route_name <- colnames(route_cnt)   #ページ名
route_M <- matrix(as.numeric(route_cnt), nrow=nrow(route_cnt), ncol(route_cnt))   #行列に変換
colnames(route_M) <- route_name


#kpi関連ページのダミー
answer <- c("answer", "survey_info")
teikei <- c("student_free", "tenshoku")
pay <- c("payment_info")

answer.v <- as.vector(rowSums(route_M[, answer]))
teikei.v <- as.vector(rowSums(route_M[, teikei]))
pay.v <- as.vector(route_M[, pay])
regist.v <- as.vector(route_M[, "regist_complete"])
pre.v <- as.vector(rowSums(route_M[, c(answer, teikei, pay)]))

#アクセス数が1000未満のアクセスのページはその他に分類
kpi.gv <- setdiff(route_name, c(answer, teikei, pay))   #不要変数を除去
route_gv <- as.matrix(route_M[, kpi.gv])
gv_name <- colnames(route_gv)

#アクセス数が1000件以上のページを抽出
index.route <- subset(1:length(gv_name), colSums(route_gv) > 1000)   

#アクセス数が1000件以上のページのダミー
route1_v <- as.matrix(route_gv[, index.route])
route1_name <- colnames(route1_v)

#その他のページのダミー
route2_v <- as.vector(rowSums(route_gv[, -index.route]))

#データを結合
page_freq1 <- data.frame(answer=answer.v, teikei=teikei.v, pay=pay.v, regist=regist.v, route1_v, other=route2_v)
page_freq2 <- data.frame(pre.v, regist=regist.v, route1_v, other=route2_v)


####閲覧企業のクレンジング####
##企業情報を結合
table(as.character(access_data_s$company_id1))
table(as.character(access_data_s$company_id2))

#企業アクセス情報をデータフレームに変換
company1 <- data.frame(no=1:nrow(access_data_s), master_id=access_data_s$company_id1)
company2 <- data.frame(no=1:nrow(access_data_s), act_id__c=access_data_s$company_id2)

#idを文字列に変換
company_vcs$act_id__c <- as.character(company_vcs$act_id__c)
company_vcs$master_id <- as.character(company_vcs$master_id)
company_all$master_id <- as.character(company_all$master_id)
company1$master_id <- as.character(company1$master_id)
company2$act_id__c <- as.character(company2$act_id__c)


#company_vcsのact_idがnullの行を削除
index.null <- subset(1:nrow(company_vcs), company_vcs$act_id__c=="")   

##企業情報と閲覧企業を結合
company1_join <- left_join(company1, company_all, by="master_id")
company2_join <- left_join(company2, company_vcs[-index.null, ], by="act_id__c")

#noで重複した行を削除する
company1_join.dp <- company1_join[!duplicated(company1_join$no), ]
company2_join.dp <- company2_join[!duplicated(company2_join$no), ]

#結合したデータを確認
company_j <- data.frame(id1=company1$master_id, id2=company2$act_id__c, route, mid1=company1_join.dp$master_id, 
                        mid2=company2_join.dp$master_id, aid=company2_join.dp$act_id__c, 
                        c1=company1_join.dp$company_name, c2=company2_join.dp$company_name)

#2つのjoinしたcompanyを1つの行列にまとめる
index.join1 <- subset(1:nrow(company1_join.dp), company1_join.dp$master_id!="")
index.join2 <- subset(1:nrow(company2_join.dp), company2_join.dp$master_id!="")

#companyを統合
company_join <- company1_join.dp
company_join[index.join2, ] <- company2_join.dp[index.join2, -2]

company_master <- data.frame(id1=company1$master_id, id2=company2$act_id__c, route, company_join)

##口コミ数ごとにダミー変数化
vcs <- paste(company_join$public_vcs_flag, company_join$public_vcs_30_flag)
v_name <- c("その他", "30件以上", "30件未満", "なし")
vcs_d <- matrix(as.numeric(table(1:length(vcs), vcs)), nrow=length(vcs), ncol=length(v_name))
colnames(vcs_d) <- v_name
vcs_d <- vcs_d[, -1]
cbind(ID, vcs_d, Z1, Z2)


##閲覧の仕方を取得
company_freq <- data.frame(id=user_id[Z2!=1], freq=company_join$master_id[Z2!=1])
company_table <- table(company_freq$id, company_freq$freq)
company_name <- colnames(company_table)
company_name <- company_name[-1]


#企業閲覧履歴を行列化
Cfreq_M0 <- matrix(as.numeric(company_table), nrow=nrow(company_table), ncol=ncol(company_table))
Cfreq_M <- Cfreq_M0[, -1]


#企業閲覧履歴から閲覧数、閲覧企業数、最大閲覧数を取得
company_sum <- rowSums(Cfreq_M)
company_max <- apply(Cfreq_M, 1, max)
company_cnt <- apply(Cfreq_M, 1, function(x) sum(x > 0))


#データを確認
#GLMのテスト
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

##セカンドセッションの指示変数
regist_time 
access_time


####コンテンツ閲覧履歴の処理####
q_no <- as.character(access_data_s$q_no)
q_no[is.na(q_no)] <- 0   #naを0にしておく

#ダミー変数化
q_no_table <- table(1:length(q_no), q_no)[, -1]
q_no_name <- colnames(q_no_table)
contents <- matrix(as.numeric(q_no_table), nrow=length(q_no), ncol=length(q_no_name))   #行列変換

#データを確認
colSums(contents)
q_no_name


##GLMテスト2
X_fun <- as.matrix(data.frame(user_id, contents) %>% 
                     dplyr::group_by(user_id) %>%
                     dplyr::summarize_each(dplyr::funs(sum)))

DATA2 <- data.frame(Y=Y_TEST, sum=company_sum, max=company_max, cnt=company_cnt, c=X_fun[, -1])
res <- glm(Y ~ ., data=DATA2[, -3], family=binomial)
summary(res)


####連続変数の処理####
#累積セッション数
session_cnt <- session_id

#登録時の登録時間
toroku <- log(access_data_s$toroku_time)
toroku[is.infinite(toroku)] <- 0


####二値変数をダミー変数に変換####
keyword <- ifelse(access_data_s$keyword=="", 0, 1)   #キーワード検索したかどうか  
sex <- ifelse(access_data_s$sex=="男性", 1, 0)   #性別
device <- ifelse(access_data_s$user_device=="sp", 1, 0)   #ユーザーデバイス
career <- ifelse(access_data_s$career=="正規社員", 1, 0)   #雇用形態
mail_follow <- ifelse(access_data_s$mail_follow==1, 1, 0)   #メールのフォロー有無
j_change <- ifelse(is.na(access_data_s$job_change_timing), 0, 1)   #転職時期が決まっているかどうか


##多値変数をダミー変数に変換
#アドミッションタイプ
admission <- matrix(as.numeric(table(1:length(user_id), access_data_s$admission_type)), nrow=length(user_id), ncol=3)


##都道府県をグループ化して、ダミー変数化
pref <- access_data_s$pref
tokyo <- c("東京都")
kanto <- c("埼玉県", "神奈川県", "千葉県", "栃木県", "茨城県", "群馬県")
kinki <- c("大阪府", "京都府", "兵庫県", "奈良県", "滋賀県", "和歌山県")
tyukyo <- c("愛知県", "三重県", "岐阜県")

#都道府県を地域ごとにグループ化
pref_o <- ifelse(pref %in% tokyo, "東京都",  ifelse(pref %in% kanto, "関東圏", ifelse(pref %in% kinki, "近畿圏",
                                                                                ifelse(pref %in% tyukyo, "中京圏", "その他"))))

#都道府県をグループ別にダミー変数化
pref_table <- table(1:length(pref), pref_o)
pref_name <- colnames(pref_table)
pref_gr <- matrix(as.numeric(pref_table), nrow=length(pref), ncol=length(pref_name))   #ダミー変数化
colnames(pref_gr) <- pref_name

colSums(pref_gr)   #データを確認


##年齢を年代に変換してダミー変数化
age <- access_data_s$age 
age_class_o <- ifelse(age < 20, "10歳代", ifelse(age >= 20 & age < 30, "20歳代", ifelse(age >= 30 & age < 40, "30歳代",
                                                                                    ifelse(age >= 40 & age < 50, "40歳代", ifelse(age >= 50 & age < 60, "50歳代", "60歳代")))))
age_table <- table(1:length(age), age_class_o)
age_name <- colnames(age_table)

#ダミー変数化
age_class <- matrix(as.numeric(age_table), nrow=length(age), ncol=length(age_name))
colnames(age_class) <- age_name


##所属業界をダミー変数化
industry_o <- access_data_s$industry
industry_table <- table(1:length(industry_o), industry_o)
industry_name <- colnames(industry_table)

industry.m <- matrix(as.numeric(industry_table), nrow=length(industry_o), ncol=length(industry_name))   #ダミー変数化
industry <- industry.m[, -1]
industry_name[-1]
colSums(industry)   #データを確認


##職業をダミー変数化
occupation_o <- access_data_s$occupation_l
occupation_table <- table(1:length(occupation_o), occupation_o)
occupation_name <- colnames(occupation_table)

occupation.m <- matrix(as.numeric(occupation_table), nrow=length(occupation_o), ncol=length(occupation_name))   #ダミー変数化
occupation <- occupation.m[, -1]
occupation_name[-1]
colSums(occupation)   #データを確認


##年収クラス
income_o <- access_data_s$income_class
income_table <- table(1:length(income_o), income_o)
income_table <- income_table[, -1]
income_name <- colnames(income_table)

income <- matrix(as.numeric(income_table), nrow=length(income_o), ncol=length(income_name))   #ダミー変数化
colSums(income)   #データの確認

##デバイス利用履歴
device_use <- matrix(0, nrow=n, ncol=3)


#ユーザーごとに利用デバイスを取得
for(i in 1:length(user_unique)){
  print(i)
  if(max(device[user_id==user_unique[i]])==1 & length(unique(device[user_id==user_unique[i]]))==1) {
    print("スマホ")
    device_use[i, 1] <- 1
    
  } else if(max(device[user_id==user_unique[i]])==0 & length(unique(device[user_id==user_unique[i]]))==1) {
    print("PC")
    device_use[i, 2] <- 1
    
  } else {
    print("両方")
    device_use[i, 3] <- 1
  }
}
colnames(device_use) <- c("sp", "PC", "両方")   #名前をつける
colSums(device_use)


####統計モデルのためのデータを作成####
####潜在変数モデリングのための説明変数を作成####
##デモグラフィック変数
Demo <- data.frame(sex, age_class, pref_gr, career, job=occupation, income=income, industry=industry, toroku=toroku,
                   ad=admission, mail=mail_follow, purpose=purpose_d, j_change=j_change)[Z2!=1, ]
demo_name <- colnames(Demo)


#ユーザーごとにデモグラフィック変数を1つにまとめる
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


#反応変数を設定
Y.demo <- data.frame(Y, ID) %>%
  dplyr::group_by(user_id) %>%
  dplyr::summarize(Y=max(Y))
Yi <- as.matrix(Y.demo)[, 2]

####生存モデルのためのパネルデータを作成####
##ユーザー登録前とあとを峻別する
regist_flag <- c()
for(i in 1:n){
  print(i)
  
  #user_idごとに登録前後のフラッグを建てる
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
  #ベクトルを結合
  regist_flag <- c(regist_flag, y_regist)
}

##ユーザー登録から24時間以内とそれ以上に分ける
diff <- access_time - regist_time   #登録時間とアクセス時間の差分
day_s <- 60*60*24   #1日を秒換算
SS_flag <- ifelse(diff >= day_s, 1, 0)


##ユーザー登録前の行動を記録
#ファーストランディンからユーザー登録までの経過日数
#ユーザー登録前にcv関連ページを踏んでいるかどうか
keika_time <- c()
keika_ind <- c()
kpi_cnt <- c()
kpi_pv <- ifelse(route %in% c(answer, teikei, pay), 1, 0)   #全体アクセスのkpiアクセス有無

for(i in 1:n){
  print(i)
  
  #ユーザーid別にファーストランディングから登録までの経過時間を記録
  index.ind <- subset(1:length(access_time), ID[, 2] == user_unique[i])
  keika <- rep(0, length(index.ind))
  p_time <- abs(min(access_time[user_id==user_unique[i]] - regist_time[user_id==user_unique[i]])/60/60/24)
  keika[1:length(keika)] <- p_time
  keika_ind <- c(keika_ind, p_time)
  keika_time <- c(keika_time, keika)   #ベクトルに結合
  
  #ユーザー登録前にcv関連ページを踏んでいるかどうかを記録
  k_cnt <- kpi_pv[user_id==user_unique[i] & regist_flag==0]
  kpi_cnt <- c(kpi_cnt, sum(k_cnt))   #ベクトルに結合
}

##ユーザー登録前の総pv数を記録
unregist_pv <- table(user_id, regist_flag)
ur_pv <- matrix(as.numeric(table(user_id, regist_flag)), nrow=length(user_unique), ncol=2)[, 1]


##ユーザー登録前のコンテンツの閲覧数
q_cnt <- table(user_id[regist_flag==0], q_no[regist_flag==0])
index_q <- subset(1:n, 1:n %in% as.numeric(rownames(q_cnt)) == FALSE)

#欠損しているidに0ベクトルを代入
cont_q <- matrix(0, nrow=n, ncol=ncol(q_cnt))

cont_q[-index_q, ] <- q_cnt
cont_q[index_q, ] <- rep(0, ncol(q_cnt))
cont_m <- cont_q[, -1]

##ユーザー登録前のページ閲覧数
pagename <- colnames(page_freq1)
page_unit <- as.numeric(as.matrix(page_freq1) %*% 1:length(pagename))

#ユーザー登録前のページの閲覧数
p_cnt <- table(c(user_id[regist_flag==0], rep(3001, ncol(page_freq1))), c(page_unit[regist_flag==0], 1:length(pagename)))
p_cnt <- p_cnt[-nrow(p_cnt), ]
index_p <- subset(1:n, 1:n %in% as.numeric(rownames(p_cnt)) == FALSE)

#欠損しているidに0ベクトルを代入
page_q <- matrix(0, nrow=n, ncol(page_freq1))

page_q[-index_p, ] <- p_cnt
page_q[index_p, ] <- rep(0, ncol(page_freq1))
colnames(page_q) <- pagename


##セッションidごとにデータをまとめてパネルデータにする
#ユーザーID、セッションごとにインデックスを記録する
session.list <- list()

#id、セッションごとにリスト化
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


#セッションIDを1行にまとめたIDを作成
ID_D <- as.matrix(ID %>% 
                    dplyr::group_by(user_id) %>%
                    dplyr::count(session_id))[, -3]

ID.v <- data.frame(no=1:nrow(ID_D), ID_D)


##セッションIDごとにアクセス日、アクセス週、アクセス時間などセッションで一意な変数を取得
#抽出する対象のデータ
TIME_DATA <- data.frame(ID, Y, Z1, Z2, regi=regist_flag, access_date, weekday, time_d, SS=SS_flag, 
                        domain=domain_d, q=contents, p=page_freq1, s_cnt=session_cnt-1, device)

#ID、セッションごとに説明変数を記録
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

#リスト形式を行列形式に変更
X.panel1 <- do.call(rbind, XZ)
C.panel1 <- do.call(rbind, Cont)
P.panel1 <- do.call(rbind, Page)


##目的、企業閲覧履歴のデータを紐付ける
#企業閲覧履歴を紐付ける
com <- matrix(0, nrow=nrow(X.panel1), ncol=3)
company_cm <- cbind(company_cnt, company_max, company_sum)

for(i in 1:n){
  index.co <- subset(1:nrow(com), ID.v$user_id==i)
  com[index.co, ] <- matrix(company_cm[i, ], nrow=length(index.co), ncol=3, byrow=TRUE)
}


#目的を紐付ける
pur <- matrix(0, nrow=nrow(X.panel1), ncol=ncol(purpose_d))

for(i in 1:n){
  index.pur <- subset(1:nrow(pur), ID.v$user_id==i)
  if(length(index.pur) > 1) {
    pur[index.pur, ] <-  matrix(purpose_d[ID$user_id==i, ][1, ], nrow=length(index.pur), ncol=ncol(purpose_d), byrow=TRUE)
  } else {
    pur[index.pur, ] <-  matrix(purpose_d[ID$user_id==i, ], nrow=length(index.pur), ncol=ncol(purpose_d), byrow=TRUE)
  } 
}


##ページ閲覧履歴、コンテンツ閲覧履歴を記録する
C.id1 <- cbind(X.panel1[, c("user_id", "session_id")], cn=C.panel1, cp=matrix(0, nrow(C.panel1), ncol(C.panel1)))
C.id2 <- cbind(C.id1, cf=matrix(0, nrow(C.panel1), ncol(C.panel1)))
P.id1 <- cbind(X.panel1[, c("user_id", "session_id")], n=P.panel1, p=matrix(0, nrow(P.panel1), ncol(P.panel1)))
P.id2 <- cbind(P.id1, pf=matrix(0, nrow(P.panel1), ncol(P.panel1)))

r1 <- 3:((ncol(P.id1)-2)/2+2)
r2 <- ((ncol(P.id1)-2)/2+3):ncol(P.id1)
r3 <- (ncol(P.id1)+1):ncol(P.id2)

for(i in 1:n){  
  print(i)
  #過去の閲覧履歴を代入
  cont_hist <- cont_m[i, ]
  page_hist <- page_q[i, ]
  index.user <- subset(1:nrow(C.id1), C.id1$user_id==i)
  user_session <- length(index.user)
  
  for(j in 1:user_session){
    C.id1[index.user[j], 13:ncol(C.id1)] <- cont_hist
    C.id2[index.user[j], 13:ncol(C.id1)] <- cont_hist
    P.id1[index.user[j], r2] <- page_hist
    P.id2[index.user[j], r2] <- page_hist
    
    #過去の閲覧履歴と現在の閲覧履歴を取得
    c.past <- C.id1[index.user[j], 13:ncol(C.id1)] 
    c.now <- C.id1[index.user[j], 3:12]
    p.past <- P.id1[index.user[j], r2] 
    p.now <- P.id1[index.user[j], r1]
    
    #初めて見るなら頻度を数える
    c.first_cnt <- abs(c.past - c.now)
    c.first_flag <- c.past-c.now < 0 & c.past==0
    c.first_cnt[c.first_flag==FALSE]  <- 0 
    
    p.first_cnt <- abs(p.past - p.now)
    p.first_flag <- p.past-p.now < 0 & p.past==0
    p.first_cnt[p.first_flag==FALSE]  <- 0
    
    #初回閲覧データを代入
    C.id2[index.user[j], (ncol(C.id1)+1):ncol(C.id2)] <- c.first_cnt
    P.id2[index.user[j], r3] <- p.first_cnt
    
    #閲覧履歴を更新
    if(j==user_session) {
      break
    } else {
      cont_hist <- C.id2[index.user[j], 13:ncol(C.id1)] + C.id2[index.user[j], 3:12]
      page_hist <- P.id2[index.user[j], r2] + P.id2[index.user[j], r1]
    } 
  }
}


##変数のスケール化
cscale <- scale(C.id2[, 13:22]) 
cp <- cscale + matrix(abs(apply(cscale, 2, min)), nrow=nrow(cscale), ncol=ncol(cscale), byrow=TRUE)
c.ones <- ifelse(cp > 0, 1, 0)   #過去に閲覧しているかどうか

pscale <- scale(P.id2[, r2])
pp <- pscale + matrix(abs(apply(pscale, 2, min)), nrow=nrow(pscale), ncol=ncol(pscale), byrow=TRUE)
p.ones <- ifelse(P.id2[, r2] > 0, 1, 0)

p_cnt.scale <- scale(rowSums(P.id2[, 33:62]))
p_cnt <- p_cnt.scale + abs(min(p_cnt.scale))

com.rate <- com[, 1:2]/com[, 3]
com.rate[is.nan(com.rate)] <- 0
com_cnt <- scale(com[, 3]) + abs(min(scale(com[, 3])))


####生存モデルとクロスセクションモデルを当てはめる####
##パネルデータを作成して生存モデルを当てはめる
#データの設定
ID.v1 <- X.panel1[, 2:3]   #IDの設定

fit_data1 <- data.frame(Y=X.panel1$Z1, X.panel1[, c(3, 9:21)], C.id2[, 3:12], C.id2[, 23:32],   #妥当 
                        P.id2[, c(69, 80, 81, 87)], P.id2[, c(9, 20, 21, 27)], pur=pur, cnt=p_cnt)

fit_data2 <- data.frame(Y=X.panel1$Z1, X.panel1[, c(3, 9:21)], C.id2[, 23:32], ones=c.ones,
                        P.id2[, c(69, 80, 81, 87)], P.id2[, c(9, 20, 21, 27)], pur=pur, cnt=p_cnt)

fit_data3 <- data.frame(Y=X.panel1$Z1, X.panel1[, c(3, 9:21)], C.id2[, 3:12], C.id2[, 23:32], ones=c.ones,
                        P.id2[, c(69, 80, 81, 87)], P.id2[, c(9, 20, 21, 27)], pur=pur, cnt=p_cnt)


ID.v2 <- ID.v1[X.panel1$Z2==0, ]
fit_data01 <- fit_data1[X.panel1$Z2==0, -c(5, 15)]
fit_data02 <- fit_data2[X.panel1$Z2==0, -c(5, 15)]
fit_data03 <- fit_data3[X.panel1$Z2==0, -c(5, 15)]

#GLMを当てはめる
fit1 <- glm(Y ~ ., data=fit_data01, family=binomial)
fit2 <- glm(Y ~ ., data=fit_data02, family=binomial)
fit3 <- glm(Y ~ ., data=fit_data03, family=binomial)

summary(fit1)
summary(fit2)
summary(fit3)


##デモグラフィック変数でGLMを当てはめる(予備分析)
#データの設定
Z.demo_add <- data.frame(Z.demo, cnt=company_cnt, sum=company_sum, device_use)

#冗長な変数を除去
index.z <- setdiff(colnames(Z.demo_add), c("X10歳代", "その他", "job.12", "ad.2", "income.6", "industry.9", "PC"))   
Zi.demo <- Z.demo_add[, index.z]

res <- glm(Y ~ ., data = data.frame(Y=Yi, Zi.demo[, -c(12:25)]), family=binomial)
summary(res)


####EMアルゴリズムでSplit_Hazard_modelを推定####
####EMアルゴリズムで利用する関数の設定####
##ロジットモデルの対数尤度の定義
logit_LL <- function(x, Y, X){
  #パラメータの設定
  b0 <- x[1]
  b1 <- x[2:(ncol(X)+1)]
  
  #効用関数の定義
  U <- b0 + as.matrix(X) %*% b1
  
  #対数尤度を計算
  Pr <- exp(U)/(1 + exp(U))   #確率の計算
  LLi <- Y*log(Pr) + (1-Y)*log(1-Pr)
  LL <- sum(LLi)
  return(LL)
}

##完全データでのSplit_hazard_modeの対数尤度
#パラメータの設定
cll <- function(b, Y, Yn, X, zpt, zk){
  #パラメータの設定
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #効用関数を定義
  U <- beta0 + as.matrix(X) %*% beta1
  
  #ロジスティックモデルの確率と対数尤度の計算
  Pr <- exp(U)/(1 + exp(U))   #確率の計算
  LLc <- Y*log(Pr) + (1-Y)*log(1-Pr)
  
  #ゼロ尤度の計算
  LLzero <- dbinom((1-Yn), 1, 1)   #ゼロ尤度の計算
  
  LL <- sum(zpt[, 1]*LLc)
  return(LL)
}


##観測データでの尤度と潜在変数zの計算
ollz <- function(b, Y, Yn, X, r, id, n, zk){
  
  #パラメータの設定
  beta0 <- b[1]
  beta1 <- b[2:(ncol(X)+1)] 
  
  #セグメントごとの効用関数を定義
  U <- beta0 + as.matrix(X) %*% beta1
  
  #ロジスティックモデルの確率と対数尤度の計算
  Pr <- exp(U)/(1 + exp(U))   #確率の計算
  LLi <- Pr^Y * (1-Pr)^(1-Y)   #尤度の計算
  LLzero <- dbinom((1-Yn), 1, 1)   #ゼロ尤度の計算
  LCo <- cbind(LLi, LLzero)   #尤度の結合
  
  #ID別に尤度の積を取る
  LLho <- matrix(0, nrow=n, ncol=zk)
  for(i in 1:n){
    if(sum(id==i)==1){
      LLho[i, ] <- LCo[id==i, ]
    } else {
      LLho[i, ] <- apply(LCo[id==i, ], 2, prod) 
    }
  }
  
  #観測データでの対数尤度
  LLo <- sum(log(apply(r * LLho, 1, sum)))
  
  #潜在変数zの計算
  z0 <- r * LLho   #潜在変数zの分子
  z1 <- z0 / rowSums(z0)   #潜在変数zの計算
  
  rval <- list(LLo=LLo, z1=z1)
  return(rval)
}

####EMアルゴリズムの設定と準備####
##EMアルゴリズムの初期値の設定
#データの設定
ID.data <- data.frame(no=1:nrow(ID.v2), ID.v2)   #IDとセッション番号
id <- ID.data$user_id
Y1 <- fit_data01[, 1]   #生存モデルの応答変数
XM <- fit_data01[, -c(1:3)]   #生存モデルの説明変数
XM[, 11:20] <- XM[, 11:20] - XM[, 21:30]
XM[, 35:38] <- XM[, 35:38] - XM[, 31:34]

Zi <- Zi.demo   #潜在モデルの説明変数
Y2 <- Yi
zpt <- matrix(0, nrow=nrow(XM), ncol=2)

#欠損処理
#Zi1 <- Zi
Zi[is.na(Zi)] <- 0
XM[is.na(XM)] <- 0
#Zi <- Zi[, -(34:42)]

#ゼロ尤度のためのyの設定
index.id <- subset(ID.data$user_id, Y1==1)
Y.zeros <- ifelse(ID.data$user_id %in% index.id, 1, 0)


####EMアルゴリズムでスプリットハザードモデルを推定####
##最良値が出るまで反復を繰り返す
rt <- 40   #繰り返し数
LL.process <- c()
LL.list <- c()
beta.list <- matrix(0, nrow=rt, ncol=ncol(XM)+1)
theta.list <- matrix(0, nrow=rt, ncol=ncol(Zi)+1)
z.list <- array(0, dim=c(n, 2, rt))

for(rp in 1:rt){
  print(rp)
  
  ##EMアルゴリズムの初期値を設定
  #生存モデルの初期値を設定
  res.surv <- glm(Y1 ~ ., data=XM, family=binomial(link="logit"))
  summary(res.surv)
  beta <- coef(res.surv)
  
  #潜在モデルの初期値を設定
  res.latent <- glm(Y2 ~ ., data=Zi, family=binomial(link="logit"))
  summary(res.latent)
  theta <- coef(res.latent)
  
  #確率の計算
  Z <- as.matrix(cbind(1, Zi)) %*% theta   #ロジットの計算
  Pr.z <- exp(Z)/(1+exp(Z))   #確率の計算
  r <- cbind(Pr.z, 1-Pr.z)   #混合率の計算
  
  ##1回目のステップのみMステップを初期値に依存しないSANNで推定
  #潜在変数zと観測データの対数尤度の初期値の設定
  
  oll <- ollz(b=beta, Y=Y1, Yn=Y.zeros, X=as.matrix(XM), r=r, id=ID.data$user_id, n=n, zk=2)
  z <- oll$z1
  LL1 <- oll$LLo
  
  #潜在変数zをパネル形式に変更
  for(i in 1:n){
    index.id <- subset(1:length(id), id==i)
    zpt[index.id, ] <- matrix(z[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
  }
  
  
  #完全データでの生存モデルの最尤推定(Mステップ)
  beta.first <- beta   #初期値を設定
  
  for(i in 1:10000){
    beta1 <- beta.first + runif(length(beta), -0.7, 0.7)
    res <- try(optim(beta1, cll, Y=Y1, Yn=Y.zeros, X=XM, zpt=zpt, zk=2, method="Nelder-Mead", 
                     hessian=FALSE, control=list(fnscale=-1, trace=TRUE, maxit=3000)), silent=TRUE)
    if(class(res) == "try-error") {next} else {break}   #エラー処理
  }
  beta <- res$par
  
  #混合率rの更新
  res.latent <- glm(z ~ ., data=Zi, family=binomial(link="logit"))
  theta <- coef(res.latent)    #潜在変数のパラメータ
  
  #確率の計算
  Z <- as.matrix(cbind(1, Zi)) %*% theta   #ロジットの計算
  Pr.z <- exp(Z)/(1+exp(Z))   #確率の計算
  r <- cbind(Pr.z, 1-Pr.z)   #混合率の計算
  
  #潜在変数zと観測データの対数尤度の初期値の設定
  oll <- ollz(b=beta, Y=Y1, Yn=Y.zeros, X=as.matrix(XM), r=r, id=ID.data$user_id, n=n, zk=2)
  z <- oll$z1
  LL1 <- oll$LLo
  
  
  #EMアルゴリズムの設定
  iter <- 0
  dl <- 100   #EMステップでの対数尤度の差の初期値を設定
  tol <- 1   
  
  ##EMアルゴリズムでスプリットハザードモデルを推定
  while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
    #潜在変数zをパネル形式に変更
    for(i in 1:n){
      index.id <- subset(1:length(id), id==i)
      zpt[index.id, ] <- matrix(z[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
    }
    
    #完全データでの生存モデルの最尤推定(Mステップ)
    res <- optim(beta, cll, Y=Y1, Yn=Y.zeros, X=XM, zpt=zpt, zk=2, method="Nelder-Mead", 
                 hessian=FALSE, control=list(fnscale=-1, maxit=2000))
    beta <- res$par
    
    #混合率rの更新
    res.latent <- glm(z ~ ., data=Zi, family=binomial(link="logit"))
    theta <- coef(res.latent)    #潜在変数のパラメータ
    Z <- as.matrix(cbind(1, Zi)) %*% theta   #ロジットの計算
    Pr.z <- exp(Z)/(1+exp(Z))   #確率の計算
    r <- cbind(Pr.z, 1-Pr.z)   #混合率の計算
    
    #潜在変数zと観測データの対数尤度を更新
    oll <- ollz(b=beta, Y=Y1, Yn=Y.zeros, X=as.matrix(XM), r=r, id=ID.data$user_id, n=n, zk=2)
    z <- oll$z1
    LL <- oll$LLo
    
    #EMアルゴリズムのパラメータの更新
    iter <- iter+1
    dl <- LL-LL1
    LL1 <- LL
    LL.process <- c(LL.process, LL)
    print(LL)
    if(abs(dl) < 5 & LL < max(LL.process)) {break}
  }
  
  #パラメータを代入
  LL.list <- c(LL.list, LL) 
  beta.list[rp, ] <- beta
  theta.list[rp, ] <- theta
  z.list[, , rp] <- z 
}

beta <- res$par

##最適なパラメータを決定
LL.list
opt.LL <- which.max(LL.list)
LL.list[opt.LL]   #提案モデルの最大対数尤度
logLik(res.surv)   #ロジスティック回帰の対数尤度

#推定されたパラメータ
beta.best <- beta.list[opt.LL, ]
theta.best <- theta.list[opt.LL, ]
z.best <- z.list[, , opt.LL]


#回帰係数に名前をつける
beta.name <- c("切片", colnames(XM))
theta.name <- c("切片", colnames(Zi))
names(beta.best) <- beta.name
names(theta.best) <- theta.name

#回帰係数を確認
round(beta.best, 3)
round(coef(res.surv), 3)
round(theta.best, 3)

occupation_name

colSums(Zi)

round(zzz <- data.frame(pre=z.best[, 1], non_pre=z.best[, 2], Y2), 3)   #潜在変数と反応変数の比較
colSums(z.best)/n   #混合率

colSums(ifelse(X[, 11:30] > 0, 1, 0))
summary(XM)

##適合度とt検定
round(res.surv$aic, 2)
round(AIC <- -2*LL.list[opt.LL] + 2*length(beta.best), 2)   ##AIC


res.latent <- glm(z.best ~ ., data=Zi, family=binomial(link="logit"))
summary(res.latent)


for(i in 1:n){
  index.id <- subset(1:length(id), id==i)
  zpt[index.id, ] <- matrix(z.best[i, ], nrow=length(index.id), ncol=2, byrow=TRUE)
}

#完全データでの生存モデルの最尤推定(Mステップ)
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

path_output <- paste("C:/Users/okuno/Desktop/", "カルテ", ".csv", sep="")   #保存パス
write.table(XZM, path_output, row.names=FALSE, col.names = TRUE,
            append=FALSE, sep=",", fileEncoding="Shift-JIS")   #書き出し


length(res.surv$fitted.values)
length(Y1)
nrow(XM)




