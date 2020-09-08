library(RSelenium)
library(RCurl)

# input parameters
topics <- c("erlotinib")
species <- "9606"
outpath <- "C:\\Users\\Grace Wu\\text_features\\R\\genes\\"

# initiate remote connection
# note: RSelenium requires docker container having run the following commands
# $ docker pull selenium/standalone-firefox:2.53.0
# $ docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.0
url <- "http://cbdm-01.zdv.uni-mainz.de/~jfontain/cms/?page_id=6"
remDr <- remoteDriver(remoteServerAddr = "192.168.99.100", port = 4445L)
remDr$open()
Sys.sleep(5)

#loop through topics of interest
for (topic in topics) { 
  remDr$navigate(url)
  webElem <- remDr$findElements("css", "iframe")
  remDr$switchToFrame(webElem[[1]])
  
  # input required fields into form
  webElem <- remDr$findElement(using = "name", value = "trset_text")
  webElem$sendKeysToElement(list(topic))
  webElem <- remDr$findElement(using = "name", value = "tset_text")
  webElem$sendKeysToElement(list(species))
  webElem <-remDr$findElement('xpath', '//*[contains(@value,"Rank it!")]')
  webElem$clickElement()
  webElem <-remDr$findElement('xpath', '//*[contains(@href,"_table.txt")]')
  webElem$clickElement()
  txt<-remDr$findElement(using='css selector',"body")$getElementText()
  
  # extract list of genes from txt output
  genes = fread(txt[[1]], sep = ' ', fill=TRUE)
  saveRDS(genes, paste(outpath, sprintf("%s.rds", topic)))
}

