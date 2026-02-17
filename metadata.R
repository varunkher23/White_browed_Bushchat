metadata = data.frame()

for (i in 1:1000){#nrow(form_data){
  if (!is.na(form_data$`data-meta-audit`[i])) {
    KEY = form_data$KEY[i]
    drive_download(form_data$`data-meta-audit`[i],path = "Metadata/odk_metadata.csv", overwrite = T)
    a = read.csv("Metadata/odk_metadata.csv")
    b = a %>%
      filter(grepl("Bird_details2",node))%>%
      separate(node, sep="/",c("Prefix1","Prefix2","Bird"))%>%
      select(Bird,latitude,longitude,accuracy)%>%
      mutate(BIRD_KEY = paste(KEY, Bird, sep = "_"))
  } else {
    b = data.frame()
  }
  metadata = rbind(metadata,b)
}
  
