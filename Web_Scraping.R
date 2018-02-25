#####Webスクレイピング#####
library(RCurl)
library(XML)
library(stringr)
library(plyr)

####整形されたデータのスクレイピング####
#urlの取得
url <- "http://www.elections.state.md.us/elections/2012/election_data/"
page_parse <- htmlParse(url, encoding="utf8")
links <- getHTMLLinks(url)
filenames <- links[str_detect(links, "_General.csv")]
filenames_list <- as.list(filenames)

#urlからデータをダウンロード
#ダウンロード用の関数の定義
downloadCSV <- function(filename, baseurl, folder){
  dir.create(folder, showWarnings=FALSE)
  fileurl <- str_c(baseurl, filename)
  if(!file.exists(str_c(folder, "/", filename))){
    download.file(fileurl, destfile=str_c(folder, "/", filename))
    Sys.sleep(1)
  }
}

#ダウンロードを実行
l_ply(filenames_list, downloadCSV, baseurl="http://www.elections.state.md.us/elections/2012/election_data/", 
      folder="elec12_maryland")
length(list.files("./elec12_maryland"))
list.files("./elec12_maryland")


####PDFで保存された選挙区地図####
#urlの取得
url <- "http://planning.maryland.gov/Redistricting/2010/legiDist.shtml"
links <- getHTMLLinks(url)
filenames <- links[str_detect(links, "2010maps/Leg/Districts_")]
filenames_list <- str_extract_all(filenames, "Districts.+pdf")

#PDFをダウンロードする関数を定義
downloadPDF <- function(filename, baseurl, folder, handle){
  dir.create(folder, showWarnings = FALSE)
  fileurl <- str_c(baseurl, filename)
  if(!file.exists(str_c(folder, "/", filename))){
    pdf_file <- getBinaryURL(fileurl, curl=handle)
    writeBin(pdf_file, str_c(folder, "/", filename))
    Sys.sleep(1)
  }
}

#ダウンロードを実行
handle <- getCurlHandle(useragent=str_c(R.version$platform, R.version$version.string, sep=", "),
                        httpheader=c(from="eddid@datacollection.com"))
#l_ply(filenames_list, downloadPDF, baseurl="planning.maryland.gov/PDF/Redistricting/2010maps/Leg/",
#      folder="elec12_maryland_maps", handle=handle)


####複数のページにアクセスするためにURLを操作####
#httpsのURLを操作
https_url <- RCurl::getURL("https://www.transparency.org/news/pressreleases/year/2010", encoding="utf8")   #httpsのurlを取得
baseurl <- htmlParse(https_url, encoding="utf8")   #httpsのurlを解析
xpath <- "//div[@id='Page']/strong[2]"
total_pages <- as.numeric(xpathSApply(baseurl, xpath, xmlValue))

#プレスリリースのページを取得
max_url <- (total_pages - 1) * 10
add_url <- str_c("/P", seq(10, max_url, 10))

#URLを取得する関数を定義
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

#関数を実行
url <- "https://www.transparency.org/news/pressreleases/year/2010"
urls_list <- getPageURLs(url)
urls_list[1:3]

##ページをダウンロード
#ダウンロードを実行する関数
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

#関数を実行
handle <- getCurlHandle()
l_ply(urls_list, dlPages, folder="tp_press_2010", handle=handle)
list.files("tp_index_2010")


##個々のプレスリリースへのリンクを特定
getPressURLs <- function(folder){
  pages_parsed <- lapply(str_c(folder, "/", dir(folder)), htmlParse)
  urls <- unlist(llply(pages_parsed, getHTMLLinks))
  press_urls <- urls[str_detect(urls, "https.+/pressrelease/")]
  press_urls_list <- as.list(press_urls)
  return(press_urls_list)
}

press_urls_list <- getPressURLs(folder="tp_index_2010")
length(press_urls_list)

#プレスリリースのダウンロード
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
