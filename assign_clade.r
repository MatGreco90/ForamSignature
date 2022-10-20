assign_clade_id<- function(x) {
  
  library(tidyverse)
  library(Biostrings)
  
  my_fast<-readDNAStringSet(x)
  
  signature_key<-tibble(
    clade_id=c('L1','L2A','L2B','L3','L4','L5A','L5B','L6','L7A','L7B','L8',
               'L9','L10A','L10B','L10C','L11','L12A','L12B','L13','L14','L15',
               'L16A','L16B','L17A','L17B','L17C','L17D','L17E','L17F','L17G',
               'L18','L19','L20','L21','L22','L23A','L23B','L24','L25','L26',
               'L27A','L27B','L28A','L28B','L29A','L29B','L29C','L31','L32',
               'L33','L34','L35','L36','L37','L38','L39','L40','L41','L42A',
               'L42B','L43'),
    signature_a=c('GACAGGCAAAAAG','CCAGGTCCGGACACATCAAGGATTGACAGGCAAAAG',
                  'CCAGGTCCGGACACATCAAGGATTGACAGGCAACAG','GACAGGCAAATAGCA',
                  'GACAGGCAACCCGAAGTG','GACAGGCAATAATAC','GACAGGCAATAATATGAATT',
                  'GACAGGCAATAATACAGAG','GACAGGCAATCATGTGGA',
                  'GACAGGCAATCGTGTGAAA','GACAGGTGTTTCGTAGGTT','GACAGGCATTAG',
                  'TGACTTTGCAAATATGCTAGTCCTTT','GCTATAGTTGAAATATGCTAGTCCTTT',
                  'TACTCTGTAAATATGCTAGTCCTTT','GACAGGCGAAAGAGT',
                  'GACAGGCGAAGGAGC','GACAGGCGAAGGAGTTTC','GACAGGTAACTGTAGTGTTA',
                  'GACAGGCGATTGTAGTTAA','GACAGGCGCAGAGTA','GACAGGCGCAGGATATC',
                  'GACAGGCGCAGGATTT','GACAGGCGCTAGTCCTGGC','GACAGGCGCTAGCGCTG',
                  'GACAGGCGCTAGCGCCGGCGATCGG','GACAGGCGCTAGGGCCTACC',
                  'GACAGGCGCTAGGGTTGGGACC','GACAGGCGCTAGCGTGGCGCAGGCCAGA',
                  'GACAGGCGCTAGCTGG','GACAGGCGCTCGCAT','GACAGGCGCTGGCA',
                  'CCAGGTTAGGACATACTGAGGATTGACAGGCGCTTGTAGTTT',
                  'GGGGAAACTTACCGGGTCCGGACACACTGAGGATTGACAGGCGTTTACT',
                  'GACAGGCTCAAAAGTTTTGGTGTTTC','GACAGGTAAAAGAAGCACA',
                  'GACAGGTAAAAGAAGTTCA','GACAGGTAAAAGTATTTACA','GACAGGTAAATGGA',
                  'GACAGGTGAATGTAGTTTGTGTTT','GACAGGTGAAAGTGGTTAA',
                  'GACAGGTGAATGTGGTTAA','GACAGGTGACATAT','GACAGGTGATATATAG',
                  'GACAGGTGATTGAGGTT','GACAGGTGATTGTGATT','GACAGGTGATTGTGGTT',
                  'GACAGGTGCTTGGCA','GACAGGTGTTGGAGAACA','GACAGGTGTTGTAAAACACA',
                  'GACAGGTGTTTAATAGT',
                  'GGGGAAACTTACCGGGTCCGGACACACTGAGGATTGACAGGTGTTTACT',
                  'GACAGGTGTTTCGAAAGGTCG','GACAGGTGTTTCGAAGGTAACT',
                  'GACAGGTTAAATAGAA','GACAGGTTCAAAGT','GACAGGTTCTGTAACTGGG',
                  'GACAGGTTCTGTGAGT',
                  'GGGGAAACTTATCAGGTCCAAACACGCTGAGGATTGACAGGTTCTTATAAGA',
                  'GGGAAATCTTACCATGGTCCAAACACGCTGAGGATTGACAGGTTCTTATAAGA',
                  'GACAGGTTTATATA')) %>% 
    group_split(clade_id) 
  
  
  find_signature= function(DF) {
    provab<-vcountPattern(DF$signature_a, my_fast,
                          max.mismatch=0, min.mismatch=0,
                          with.indels=FALSE, fixed=TRUE,
                          algorithm="auto")
   pos_signature<-which(provab %in% 1)
     if (length(pos_signature)!=0){
       fasta_subset<-my_fast[pos_signature]
       
       dnadf<-as.data.frame(fasta_subset)
      
       dnadf$seq_name<-rownames(dnadf)
       colnames(dnadf)<-c('sequence','ASV_id')
       
       dnadf$assigned_clade<-DF$clade_id
       final_df<-dnadf[,c(2,1,3)]
        write.csv(final_df, paste0('assigned_',unique(DF$clade_id),".csv"),row.names=FALSE)
     }
   }
  lapply(signature_key, find_signature)
}