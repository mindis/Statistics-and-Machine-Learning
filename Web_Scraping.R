#####Web�X�N���C�s���O#####
library(RCurl)
library(XML)
library(stringr)
library(plyr)

####���`���ꂽ�f�[�^�̃X�N���C�s���O####
#url�̎擾
url <- "http://www.elections.state.md.us/elections/2012/election_data/"
page_parse <- htmlParse(url, encoding="utf8")
links <- getHTMLLinks(url)
filenames <- links[str_detect(links, "_General.csv")]
filenames_list <- as.list(filenames)

#url����f�[�^���_�E�����[�h
#�_�E�����[�h�p�̊֐��̒�`
downloadCSV <- function(filename, baseurl, folder){
  dir.create(folder, showWarnings=FALSE)
  fileurl <- str_c(baseurl, filename)
  if(!file.exists(str_c(folder, "/", filename))){
    download.file(fileurl, destfile=str_c(folder, "/", filename))
    Sys.sleep(1)
  }
}

#�_�E�����[�h�����s
l_ply(filenames_list, downloadCSV, baseurl="http://www.elections.state.md.us/elections/2012/election_data/", 
      folder="elec12_maryland")
length(list.files("./elec12_maryland"))
list.files("./elec12_maryland")


####PDF�ŕۑ����ꂽ�I����n�}####
#url�̎擾
url <- "http://planning.maryland.gov/Redistricting/2010/legiDist.shtml"
links <- getHTMLLinks(url)
filenames <- links[str_detect(links, "2010maps/Leg/Districts_")]
filenames_list <- str_extract_all(filenames, "Districts.+pdf")

#PDF���_�E�����[�h����֐����`
downloadPDF <- function(filename, baseurl, folder, handle){
  dir.create(folder, showWarnings = FALSE)
  fileurl <- str_c(baseurl, filename)
  if(!file.exists(str_c(folder, "/", filename))){
    pdf_file <- getBinaryURL(fileurl, curl=handle)
    writeBin(pdf_file, str_c(folder, "/", filename))
    Sys.sleep(1)
  }
}

#�_�E�����[�h�����s
handle <- getCurlHandle(useragent=str_c(R.version$platform, R.version$version.string, sep=", "),
                        httpheader=c(from="eddid@datacollection.com"))
#l_ply(filenames_list, downloadPDF, baseurl="planning.maryland.gov/PDF/Redistricting/2010maps/Leg/",
#      folder="elec12_maryland_maps", handle=handle)


####�����̃y�[�W�ɃA�N�Z�X���邽�߂�URL�𑀍�####
#https��URL�𑀍�
https_url <- RCurl::getURL("https://www.transparency.org/news/pressreleases/year/2010", encoding="utf8")   #https��url���擾
baseurl <- htmlParse(https_url, encoding="utf8")   #https��url�����
xpath <- "//div[@id='Page']/strong[2]"
total_pages <- as.numeric(xpathSApply(baseurl, xpath, xmlValue))

#�v���X�����[�X�̃y�[�W���擾
max_url <- (total_pages - 1) * 10
add_url <- str_c("/P", seq(10, max_url, 10))

#URL���擾����֐����`
getPageURLs <- function(url){
  https_url <- RCurl::getURL(url, encoding="utf8")
  baseurl <- htmlParse(https_url)
  total_pages <- as.numeric(xpathSApply(baseurl, "//div[@id='Page']/strong[2]", xmlValue))
  max_url <- (total_pages - 1) * 10
  add_url <- str_c("/P", seq(10, max_url, 10))
  urls_list <- as.list(str_c(url, add_url))
  urls_list[length(urls_list) + 1] <- url
  return(urls_list)
}

#�֐������s
url <- "https://www.transparency.org/news/pressreleases/year/2010"
urls_list <- getPageURLs(url)
urls_list[1:3]

##�y�[�W���_�E�����[�h
#�_�E�����[�h�����s����֐�
dlPages <- function(pageurl, folder, handle){
  dir.create(folder, showWarnings=FALSE)
  page_name <- str_c(str_extract(pageurl, "/P.+"), ".html")
  if(is.na(page_name)){page_name <- "/base.html"}
  if(!file.exists(str_c(folder, "/", page_name))){
    content <- try(getURL(pageurl, curl=handle))
    write(content, str_c(folder, "/", page_name))
    Sys.sleep(1)
  }
}

#�֐������s
handle <- getCurlHandle()
l_ply(urls_list, dlPages, folder="tp_press_2010", handle=handle)
list.files("tp_index_2010")


##�X�̃v���X�����[�X�ւ̃����N�����
getPressURLs <- function(folder){
  pages_parsed <- lapply(str_c(folder, "/", dir(folder)), htmlParse)
  urls <- unlist(llply(pages_parsed, getHTMLLinks))
  press_urls <- urls[str_detect(urls, "https.+/pressrelease/")]
  press_urls_list <- as.list(press_urls)
  return(press_urls_list)
}

press_urls_list <- getPressURLs(folder="tp_index_2010")
length(press_urls_list)

#�v���X�����[�X�̃_�E�����[�h
dlPress <- function(press_url, folder, handle){
  dir.create(folder, showWarnings=FALSE)
  press_filename <- str_c(str_extract(press_url, "[^//][[:alnum:]_.]+$"), ".html")
  if(!file.exists(str_c(folder, "/", press_filename))){
    content <- try(getURL(press_url, curl=handle))
    write(content, str_c(folder, "/", press_filename))
    Sys.sleep(1)
  }
}

handle <- getCurlHandle()
l_ply(press_urls_list, dlPress, folder="tp_press_2010", handle=handle)
length(list.files("tp_press_2010"))