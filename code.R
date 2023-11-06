#Carlos Eugenio Saldaña Tijerina

#Librerias a utilizar
library(seqinr)
library(ape)
library(Biostrings)
library(adegenet)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2)



virus <- c(  "MN985325", "MN908947", "MT012098", "MT320538", "MT394864" , "MT126808" , "LC529905", "MT039890" , "MT066156" , "LR813996" , "MT510643" , "MT327745" , "MT198651" , "MT192772" , "MT007544" , "MW553294" , "MT396266" , "MT281530" , "MT810752" , "MZ026853")

virus_sequences <- read.GenBank(virus)


alemania <- read.fasta("ALEMANIA.fasta")
cat("Tamaño de la secuencia(Alemania):", length(alemania[[1]]))

argentina <- read.fasta("ARGENTINA.fasta")
cat("Tamaño de la secuencia(Argentina):", length(argentina[[1]]))

australia <- read.fasta("AUSTRALIA.fasta")
cat("Tamaño de la secuencia(Australia):", length(australia[[1]]))

brazil <- read.fasta("BRAZIL.fasta")
cat("Tamaño de la secuencia(Brazil):", length(brazil[[1]]))

china <- read.fasta("CHINA.fasta")
cat("Tamaño de la secuencia(China):", length(china[[1]]))

espana <- read.fasta("ESPANA.fasta")
cat("Tamaño de la secuencia(España):", length(espana[[1]]))

francia <- read.fasta("FRANCIA.fasta")
cat("Tamaño de la secuencia(Francia):", length(francia[[1]]))


india <- read.fasta("INDIA.fasta")
cat("Tamaño de la secuencia(India):", length(india[[1]]))

indonesia <- read.fasta("INDONESIA.fasta")
cat("Tamaño de la secuencia(Indonesia):", length(indonesia[[1]]))


iran<- read.fasta("IRAN.fasta")
cat("Tamaño de la secuencia(Iran):", length(iran[[1]]))

italia <- read.fasta("ITALIA.fasta")
cat("Tamaño de la secuencia(Italia):", length(italia[[1]]))

japon<- read.fasta("JAPON.fasta")
cat("Tamaño de la secuencia(Japon):", length(japon[[1]]))

mexico<- read.fasta("MEXICO.fasta")
cat("Tamaño de la secuencia(Mexico):", length(mexico[[1]]))

nether<- read.fasta("NETHERLANDS.fasta")
cat("Tamaño de la secuencia(Netherlands):", length(nether[[1]]))

rusia<- read.fasta("RUSSIA.fasta")
cat("Tamaño de la secuencia(Russia):", length(rusia[[1]]))

skorea <- read.fasta("S_KOREA.fasta")
cat("Tamaño de la secuencia(South Korea):", length(skorea[[1]]))

turquia <- read.fasta("TURQUIA.fasta")
cat("Tamaño de la secuencia(Turquia):", length(turquia[[1]]))

uk <- read.fasta("UK.fasta")
cat("Tamaño de la secuencia(United Kingdom):", length(uk[[1]]))

usa<- read.fasta("USA.fasta")
cat("Tamaño de la secuencia(USA):", length(usa[[1]]))

vietnam <- read.fasta("VIETNAM.fasta")
cat("Tamaño de la secuencia(Vietnam):", length(vietnam[[1]]))


Alemania_count <- count(alemania[[1]], 1)
Argentina_count <- count(argentina[[1]],1)
Australia_count <- count(australia[[1]],1)
Brazil_count<- count(brazil[[1]],1)
China_count <- count(china[[1]],1)
Espana_count <- count(espana[[1]],1)
Francia_count <- count(francia[[1]],1)
India_count <- count(india[[1]],1)
Indonesia_count <- count(indonesia[[1]],1)
Iran_count <- count(iran[[1]],1)
Italia_count <- count(italia[[1]],1)
Japon_count <- count(japon[[1]],1)
Mexico_count <- count(mexico[[1]],1)
Nether_count <- count(nether[[1]],1)
Rusia_count <- count(rusia[[1]],1)
Skorea_count <- count(skorea[[1]],1)
Turquia_count <- count(turquia[[1]],1)
UK_count <- count(uk[[1]],1)
USA_count <- count(usa[[1]],1)
Vietnam_count <- count(vietnam[[1]],1)


counts_matrix <- rbind(Alemania_count, Argentina_count, Australia_count,Brazil_count,China_count,
                       Espana_count, Francia_count,India_count, Indonesia_count,Iran_count,
                       Italia_count, Japon_count, Mexico_count, Nether_count,
                       Rusia_count,Skorea_count,Turquia_count,UK_count, USA_count, Vietnam_count) 





barplot(counts_matrix, beside = TRUE, col = c("blue", "red", "green", "orange", "purple", "#FAEBD7", "#8B8378", 
                                              "#A52A2A", "#98F5FF", "#2F4F4F", "darkorange4", "#CDC673", 
                                              "#EEA2AD", "olivedrab4", "#CD7054", "#CDB5CD", "#8B7E66", "#FFFF00", 
                                              "#8B5A2B", "#EEE5DE"),
        
)



legend_labels <- c("Alemania", "Argentina", "Australia", "Brazil", "China", 
                           "España", "Francia", "India", "Indonesia", "Iran",
                           "Italia", "Japon", "Mexico", "Netherlands", "Rusia",
                           "Korea del Sur", "Turquia", "United Kingdom", "United States",
                           "Vietnam")
legend_colors <- c("blue", "red", "green", "orange", "purple", "#FAEBD7", "#8B8378", 
                           "#A52A2A", "#98F5FF", "#2F4F4F", "darkorange4", "#CDC673", 
                           "#EEA2AD", "olivedrab4", "#CD7054", "#CDB5CD", "#8B7E66", "#FFFF00", 
                           "#8B5A2B", "#EEE5DE")
        
par(mar=c(5, 6, 4, 2) + 0.1) # adjust the margins
plot.new() # create an empty plot
legend("center", legend = legend_labels, fill = legend_colors, ncol = 2)


write.dna(virus_sequences,  file ="virus_seqs.fasta", format = "fasta", append =
            FALSE, nbcol = 6, colsep = " ", colw = 10)

virus_seq_no_alineadas <- readDNAStringSet("virus_seqs.fasta", format = "fasta")
virus_seq_no_alineadas <- OrientNucleotides(virus_seq_no_alineadas)
virus_align_seqs <- AlignSeqs(virus_seq_no_alineadas)
writeXStringSet(virus_align_seqs, file = "virus_align_seq.fasta")
virus_aligned <- read.alignment("virus_align_seq.fasta", format = "fasta") 
matriz_distancia <- dist.alignment(virus_aligned, matrix = "similarity")
temp <- as.data.frame(as.matrix(matriz_distancia))

virus_filogenetico <- nj(matriz_distancia)
class(virus_filogenetico)
new_names = c("USA","CHN","IND","FRA","DEU","BRA","JPN","KOR","ITA","GBR","RUS","TUR","ESP","VNM","AUS","ARG","NDL","IRA","MEX","IDN")
virus_filogenetico$tip.label = new_names

virus_plot_filogenetico <- ladderize(virus_filogenetico)
ggtree(virus_plot_filogenetico) +
  geom_tiplab(size = 3) +
  theme_tree2() +
  theme(legend.position = "none", text = element_text(size = 6)) +
  ggtitle("SARS-CoV-2 en los 20 países con más casos reportados")

