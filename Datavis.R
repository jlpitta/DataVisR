#R e visualização de dados através de exemplos
#Por: Túlio Campos (Fiocruz-PE)

#Instalação de bibliotecas#

#Instala o gerenciador de pacotes do Bioconductor
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#Bibliotecas do Bioconductor
#BiocManager::install("ggtree")
#BiocManager::install("treeio")
#BiocManager::install("ggjoy")
#BiocManager::install("phyloseq")

#Bibliotecas do repositório do R (CRAN)
#install.packages("tidyverse")
#install.packages("geobr")
#install.packages("TDbook")
#install.packages("plotly")

#Diretório corrente
getwd()

#Lista arquivos no diretório
dir()

#Seta diretório de trabalho
setwd("C:/DataVis")

##############

#Carregando todas as bibliotecas que usaremos

library("ggplot2")
library("plotly")
library("geobr")
library("MASS")
library("maps")
library("ggtree")
library("TDbook")
library("treeio")
library("aplot")
library("phyloseq")
library("ggjoy")
library("dplyr")


#Revisão de R 

#Atribuição
x <- 5
x

#Checagem de tipo de dados
typeof(TRUE)
typeof(x)
class(x)
str(x)

#Lógico
typeof(TRUE)
x>5
x<5
x==5

#Lista
a_list <- list(title = "Numbers", numbers = 1:10, data = TRUE )
a_list
a_list$title
a_list[[2]]
a_list[[2]][10]

#Checagem de tipo de dados
typeof(x)
class(x)

#Ajuda
help(typeof)
help(class)
??plot

#Vetor
n <- c(1,2,3)
n

#Vetor com sequência/operação no vetor
n <- c(13:19)
n
n + 5

#Seleção de dados no vetor
n[4]
new_n <- n[c(1,4)] 
new_n

#Tamanho
length(new_w)

#Filtro
n[c(-1, -4)] 

#Erro - não existe
n[10]

#NA
is.na(n[10])

#Erro - inexistente - R é 'case sensitive'!
N 

#Conversão de tipo
as.character(n)
as.integer("1")

#if/else
if(x==5) {
  print("É igual")
} else {
  print("Não é igual")
}

#loop for
for(i in 1:5) {
  print("Teste")
  print(i)
}

#while loop
x<-1

while(x<5) {
  print(x)
  x<-x+1
}

#Vetores, Matrizes, Array (coleção de matrizes)

#Dataframes

Student_Name <- c("Amy", "Bob", "Chuck", "Daisy", "Ellie", "Frank", 
                  "George", "Helen")
Age <- c(27, 55, 34, 42, 20, 27, 34, 42)
Gender <- c("F", "M", "M", "F", "F", "M", "M", "F")
GPA <- c(3.26, 3.75, 2.98, 3.40, 2.75, 3.32, 3.68, 3.97)
nsc <- data.frame(Student_Name, Age, Gender, GPA)   
nsc 

dim(nsc)
is.na(nsc)

nsc <- nsc[!is.na(nsc)]

#Importando dados de arquivos de dados estruturados/planilhas

#read.csv("Arquivo.csv", sep=",") #Atenção para o separador
#read_excel("Arquivo.xlsx", sheet = 2)
#View(Arquivo)
#head(Arquivo)
#tail(Arquivo)

#Selecionando dados do dataframe
names(nsc)
nsc[c("Age")]   
nsc[2]
nsc[2:3]
nsc[c("Age", "Gender")]
nsc[2,] 
nsc[2:3, c("Age")] 
nsc[2:3, 2:3]
nsc$Age  
nsc$Gender  
table(nsc$Age) #counts
table(nsc$Age, nsc$Gender)  

#Operando com dataframes no dataset ChickWeight
?ChickWeight
head(ChickWeight,3)
tail(ChickWeight,5)
class(ChickWeight)
typeof(ChickWeight)

chick0 <- ChickWeight[ChickWeight$Time == 0, ]
chick0

chick0_ascend <- chick0[order(chick0$weight), ]
head(chick0_ascend, 15)

chick0_descend <- chick0[order(-chick0$weight), ]
head(chick0_descend, 15)

chick0_diet <- chick0[order(-chick0$Diet), ] 

chick0_arr <- chick0[order(chick0$Diet, chick0$weight), ]
chick0_arr

chick0_arr2 <- chick0[order(chick0$Diet, -chick0$weight), ]
chick0_arr2

names(ChickWeight)[2] <- "Days"
head(ChickWeight)

names(ChickWeight)[names(ChickWeight) == "Days"] <- "Time"
head(ChickWeight)

head(ChickWeight, 10)
ChickWeight[7, 1] = 16
head(ChickWeight, 10)

ChickWeight$weight[7] = 106
head(ChickWeight, 10)


#Tibble dataframes

data <- data.frame(a = 1:3, b = letters[1:3], c = Sys.Date() - 1:3)
data

as_tibble(data)

data %>% head(2)


#Trabalhando com gráficos

#Barplot/barras padrão do R

Browser <- c("Chrome", "Edge", "Firefox", "IE",
             "Opera", "Safari", "Others")
Users <- c(2502.4, 150.78, 395.83, 238.05, 86.49, 387.65, 134.8)
IB <- data.frame(Browser, Users)
IB

barplot(IB$Users)

IB_asc <- IB[order(IB$Users),]

png("rplot1.png", width = 350, height = 350)
barplot(IB_asc$Users, 
        main = "2018 Internet Browser Users (in millions)",
        xlab = "Internet Browser", 
        ylab = "Users", 
        names.arg = IB_asc$Browser)
dev.off()

#Usando ggplot2
ggplot(IB, aes(Browser, Users)) + 
  geom_bar(stat = "identity")

ggplot(data = IB, aes(x = Browser, y = Users)) +
  geom_bar(stat = "identity") + theme_minimal()

#theme_dark()
#theme_classic()

ggplot(data = IB, aes(x = Browser, y = Users)) + 
  geom_bar(stat = "identity") +
  coord_flip()

ggplot(data = IB, aes(x = Browser, y = Users, fill = Browser)) +
  geom_bar(stat = "identity")

#Salvando no ggplot
ggsave("rplot2.png")

#Pizza

Cities <- c("New York", "Los Angeles", "Chicago", "Houston", "Phoenix", 
            "Philadelphia", "San Antonio", "San Diego", "Dallas", "San Jose")
Population <- c(8.60, 4.06, 2.68, 2.40, 2.71, 1.58, 1.57, 1.45, 1.40, 1.03 )
top_ten <- data.frame(Cities, Population)
top_ten

pie(top_ten$Population, labels = Cities)

top_ten_descend <- top_ten[order(-top_ten$Population), ]
pie(top_ten_descend$Population, 
    labels = Cities, 
    main = "Most Populous US Cities in 2019 (in millions)")

ggplot(data = top_ten, aes(x = "", y = Population, fill = Cities)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y")

ggplot(data = top_ten, aes(x = "", y = -Population, 
                           fill = reorder(Cities, -Population))) + 
  geom_bar(stat = "identity", color = "black") + 
  labs(title = "Most Populous US Cities in 2019 (in millions)") +
  coord_polar("y") +
  theme_void()

#Histograma

help(faithful)
head(faithful)
ggplot(data = faithful, aes(x = waiting)) + 
  geom_histogram()

ggplot(data = faithful, aes(x = waiting)) + 
  geom_histogram(binwidth = 5, color = "black", fill = "white") + 
  labs(title = "IQ Score", x = "IQ Score", y = "Frequency")


#Boxplots

boxplot(ldeaths, 
        main = "Monthly Deaths from Lung Diseases in the UK", 
        ylab = "Number of Deaths")


ggplot(data = UScereal, aes(y = sodium)) + 
  geom_boxplot() + 
  scale_x_discrete() +
  labs(title = "Sodium Content in One Cup of US Cereal",
       y = "Sodium Content (in milligrams)")

boxplot(chickwts$weight ~ chickwts$feed, 
        main = "Chick Weight by Feed Type", 
        xlab = "Feed Type",
        ylab = "Weight (in grams)")

#Scatterplot

ggplot(data = Cars93, aes(x = Weight, y = MPG.city)) +
  geom_point() +
  labs(title = "Scatterplot of Weight of Car vs City MPG",
       x = "Weight of Car (in pounds)",
       y = "City Miles per Gallon")

ggplot(data = Cars93, aes(x = Weight, y = MPG.city)) +
  geom_point(color = "red") +
  geom_smooth(method = "lm", se = FALSE, col = "purple") +
  labs(title = "Scatterplot of Weight of Car vs City MPG",
       x = "Weight of Car (in pounds)",
       y = "City Miles per Gallon")

ggplot(data = Melanoma, aes(x = age, y = time, color = sex)) +
  geom_point() +
  labs(title = "Survival Time from Malignant Melanoma",
       x = "Age (in years)",
       y = "Survival Time (in days)")

#############Plotly

#Scatter

fig <- plot_ly(data = iris, x = ~Sepal.Length, y = ~Petal.Length, color = ~Species)

fig

d <- diamonds[sample(nrow(diamonds), 1000), ]

fig <- plot_ly(
  d, x = ~carat, y = ~price,
  color = ~carat, size = ~carat
)

fig

#Linha

trace_0 <- rnorm(100, mean = 5)
trace_1 <- rnorm(100, mean = 0)
trace_2 <- rnorm(100, mean = -5)
x <- c(1:100)

data <- data.frame(x, trace_0, trace_1, trace_2)

fig <- plot_ly(data, x = ~x, y = ~trace_0, name = 'trace 0', type = 'scatter', mode = 'lines') 
fig <- fig %>% add_trace(y = ~trace_1, name = 'trace 1', mode = 'lines+markers') 
fig <- fig %>% add_trace(y = ~trace_2, name = 'trace 2', mode = 'markers')

fig

#Barras e barras empilhadas

Animals <- c("giraffes", "orangutans", "monkeys")
SF_Zoo <- c(20, 14, 23)
LA_Zoo <- c(12, 18, 29)
data <- data.frame(Animals, SF_Zoo, LA_Zoo)

fig <- plot_ly(data, x = ~Animals, y = ~SF_Zoo, type = 'bar', name = 'SF Zoo')
fig <- fig %>% add_trace(y = ~LA_Zoo, name = 'LA Zoo')
fig <- fig %>% layout(yaxis = list(title = 'Count'), barmode = 'group')

fig

#

Animals <- c("giraffes", "orangutans", "monkeys")
SF_Zoo <- c(20, 14, 23)
LA_Zoo <- c(12, 18, 29)
data <- data.frame(Animals, SF_Zoo, LA_Zoo)

fig <- plot_ly(data, x = ~Animals, y = ~SF_Zoo, type = 'bar', name = 'SF Zoo')
fig <- fig %>% add_trace(y = ~LA_Zoo, name = 'LA Zoo')
fig <- fig %>% layout(yaxis = list(title = 'Count'), barmode = 'stack')

fig

#Pizza

fig <- plot_ly()
fig <- fig %>% add_pie(data = count(diamonds, cut), labels = ~cut, values = ~n,
                       name = "Cut", domain = list(x = c(0, 0.4), y = c(0.4, 1)))
fig <- fig %>% add_pie(data = count(diamonds, color), labels = ~color, values = ~n,
                       name = "Color", domain = list(x = c(0.6, 1), y = c(0.4, 1)))
fig <- fig %>% add_pie(data = count(diamonds, clarity), labels = ~clarity, values = ~n,
                       name = "Clarity", domain = list(x = c(0.25, 0.75), y = c(0, 0.6)))
fig <- fig %>% layout(title = "Pie Charts with Subplots", showlegend = F,
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
fig

#Bubble/Bolhas

data <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/school_earnings.csv")

fig <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', mode = 'markers', color = ~Gap, colors = 'Reds',
               marker = list(size = ~Gap, opacity = 0.5))
fig <- fig %>% layout(title = 'Gender Gap in Earnings per University',
                      xaxis = list(showgrid = FALSE),
                      yaxis = list(showgrid = FALSE))

fig

###Trabalhando com árvores

#Tipos de árvore

ggtree(tree)
ggtree(tree, layout="slanted") 
ggtree(tree, layout="circular")
ggtree(tree, layout="fan", open.angle=120)
ggtree(tree, layout="equal_angle")
ggtree(tree, layout="daylight")
ggtree(tree, branch.length='none')
ggtree(tree, branch.length='none', layout='circular')
ggtree(tree, layout="daylight", branch.length = 'none')

#Exemplo1

tr <- rtree(10)

d1 <- data.frame(
  label = c(tr$tip.label[sample(5, 5)], "A"),
  value = sample(1:6, 6))

d2 <- data.frame(
  label = rep(tr$tip.label, 5),
  category = rep(LETTERS[1:5], each=10),
  value = rnorm(50, 0, 3)) 

g <- ggtree(tr) + geom_tiplab(align=TRUE) + hexpand(.01)

p1 <- ggplot(d1, aes(label, value)) + geom_col(aes(fill=label)) + 
  geom_text(aes(label=label, y= value+.1)) +
  coord_flip() + theme_tree2() + theme(legend.position='none')

p2 <- ggplot(d2, aes(x=category, y=label)) + 
  geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
  theme_minimal() + xlab(NULL) + ylab(NULL)

p2 %>% insert_left(g) %>% insert_right(p1, width=.5) 

###

file <- system.file("extdata/BEAST", "beast_mcc.tree", 
                    package="treeio")
beast <- read.beast(file)
ggtree(beast) + 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=.3) + 
  geom_range("length_0.95_HPD", color='red', size=2, alpha=.5) + 
  geom_text2(aes(label=round(as.numeric(posterior), 2), 
                 subset=as.numeric(posterior)> 0.9, 
                 x=branch), vjust=0) 


#Exemplo2

tree <- read.newick(text=text_RMI_tree, node.label = "support")
root <- rootnode(tree)  
ggtree(tree, color="black", size=1.5, linetype=1,  right=TRUE) + 
  geom_tiplab(size=4.5, hjust = -0.060) +  xlim(0, 0.09) + 
  geom_point2(aes(subset=!isTip & node != root, 
                  fill=cut(support, c(0, 700, 900, 1000))), 
              shape=21, size=4) + 
  theme_tree(legend.position=c(0.2, 0.2)) + 
  scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                    name='Bootstrap Percentage(BP)', 
                    breaks=c('(900,1e+03]', '(700,900]', '(0,700]'), 
                    labels=expression(BP>=90,70 <= BP * " < 90", BP < 70))


#Exemplo3

mytree <- tree_treenwk_30.4.19

tiplab <- mytree$tip.label
cls <- tiplab[grep("^ch", tiplab)] 
labeltree <- groupOTU(mytree, cls)

p <- ggtree(labeltree, aes(color=group, linetype=group), layout="circular") +
  scale_color_manual(values = c("#efad29", "#63bbd4")) +
  geom_nodepoint(color="black", size=0.1) +
  geom_tiplab(size=2, color="black")

p2 <- flip(p, 136, 110) %>% 
  flip(141, 145) %>% 
  rotate(141) %>% 
  rotate(142) %>% 
  rotate(160) %>% 
  rotate(164) %>% 
  rotate(131)

#Agrupar V and II, colorindo
dat <- data.frame(
  node = c(110, 88, 156,136),
  fill = c("#229f8a", "#229f8a", "#229f8a", "#f9311f")
)
p3 <- p2 +
  geom_hilight(
    data = dat,
    mapping = aes(
      node = node,
      fill = I(fill)
    ),
    alpha = 0.2,
    extendto = 1.4
  )

#Adicionar label na expansão
p4 <- p3 +
  geom_cladelab(
    node = 113,
    label = "Avian-specific expansion",
    align = TRUE,
    angle = -35,
    offset.text = 0.05,
    hjust = "center",
    fontsize = 2,
    offset = .2,
    barsize = .2
  )

#Adicionar labels pros boostraps >50  
p5 <- p4 +
  geom_nodelab(
    mapping = aes(
      x = branch,
      label = label,
      subset = !is.na(as.numeric(label)) & as.numeric(label) > 50
    ),
    size = 2,
    color = "black",
    nudge_y = 0.6
  )

### Putting labels on the subgroups 
p6 <- p5 +
  geom_cladelab(
    data = data.frame(
      node = c(114, 121),
      name = c("Subgroup A", "Subgroup B")
    ),
    mapping = aes(
      node = node,
      label = name
    ),
    align = TRUE,
    offset = .05,
    offset.text = .03,
    hjust = "center",
    barsize = .2,
    fontsize = 2,
    angle = "auto",
    horizontal = FALSE
  ) +
  theme(
    legend.position = "none",
    plot.margin = grid::unit(c(-15, -15, -15, -15), "mm")
  )
print(p6)

#Exemplo4

btrees <- read.tree(system.file("extdata/RAxML", 
                                "RAxML_bootstrap.H3", 
                                package="treeio")
)
ggtree(btrees) + facet_wrap(~.id, ncol=10)

#Exemplo5

beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)
beast_tree

p1 <- ggtree(beast_tree, mrsd='2013-01-01') + theme_tree2() +
  labs(caption="Divergence time")
p2 <- ggtree(beast_tree, branch.length='rate') + theme_tree2() +
  labs(caption="Substitution rate")

mlcfile <- system.file("extdata/PAML_Codeml", "mlc", package="treeio")
mlc_tree <- read.codeml_mlc(mlcfile)
p3 <- ggtree(mlc_tree) + theme_tree2() +
  labs(caption="nucleotide substitutions per codon")
p4 <- ggtree(mlc_tree, branch.length='dN_vs_dS') + theme_tree2() +
  labs(caption="dN/dS tree")

beast_tree2 <- rescale_tree(beast_tree, branch.length='rate')
ggtree(beast_tree2) + theme_tree2()

ggdensitree(btrees, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=2) + hexpand(.35)


#Exemplo6

#Carrega `tree_nwk`, `df_info`, `df_alleles`, e `df_bar_data` do 'TDbook'

class(tree_nwk)
tree <- tree_nwk
snps <- df_alleles
snps_strainCols <- snps[1,] 
snps<-snps[-1,] # remove nomes de strains
colnames(snps) <- snps_strainCols #adiciona strains como nome de colunas

gapChar <- "?"
#Transpor dataframe
snp <- t(snps)
#Função para identificar SNP
lsnp <- apply(snp, 1, function(x) {
  x != snp[1,] & x != gapChar & snp[1,] != gapChar
})
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]

#Cria o objeto para visualização da árvore 
p <- ggtree(tree) 

# Adiciona informações das amostras
# colorindo por localidade
p <- p %<+% df_info + geom_tippoint(aes(color=location))

## Adiciona informações de SNPs e fenótipos usando gráfico de pontos/barras
## alinha com a estrutura da árvore
p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = location), shape = '|') +
  geom_facet(panel = "Trait", data = df_bar_data, geom = geom_col, 
             aes(x = dummy_bar_value, color = location, 
                 fill = location), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))

#Exemplo7

beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
colnames(genotype) <- sub("\\.$", "", colnames(genotype))
p <- ggtree(beast_tree, mrsd="2013-01-01") + 
  geom_treescale(x=2008, y=1, offset=2) + 
  geom_tiplab(size=2)

#Integrando a árvore do ggtree com gheatmap
gheatmap(p, genotype, offset=5, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
                    values=c("steelblue", "firebrick", "darkgreen"), name="genotype")

#Exemplo8

p <- ggtree(beast_tree, mrsd="2013-01-01") + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()
gheatmap(p, genotype, offset=8, width=0.6, 
         colnames=FALSE, legend_title="genotype") +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))

############### MAPAS 

#Exemplo usando geobr

datasets <- list_geobr()

datsets 

states <- read_state(
  year=2020, 
  showProgress = FALSE
)

# Remove eixos
no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

# Lê o dados de expectativa de vida e salva em dataframe
df <- utils::read.csv(system.file("extdata/br_states_lifexpect2017.csv", package = "geobr"), encoding = "UTF-8")

states$name_state <- tolower(states$name_state)
df$uf <- tolower(df$uf)

# Junta bases de dados
states <- dplyr::left_join(states, df, by = c("name_state" = "uf"))

ggplot() +
  geom_sf(data=states, aes(fill=ESPVIDA2017), color= NA, size=.15) +
  labs(subtitle="Life Expectancy at birth, Brazilian States, 2014", size=8) +
  scale_fill_distiller(palette = "Blues", name="Life Expectancy", limits = c(65,80)) +
  theme_minimal() + 
  no_axis

#Mais um exemplo com ggplot - Rio de Janeiro

idh=read.csv("https://www.dropbox.com/s/xvapy4b9ugv27zi/idh.csv?dl=1",header=TRUE, sep=";")

idh$idh=as.numeric(gsub(",", ".", gsub("\\.", "", idh$idh)))

mun_rj <- read_municipality(code_muni=33, year=2017)
all_mun_rj = left_join(mun_rj, idh, by='name_muni')

no_axis <- theme(axis.title=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank())

ggplot() +
  geom_sf(data=all_mun_rj, aes(fill=idh), color= NA, size=.15)+
  labs(title="IDH dos Municipíos do RJ",
       caption='Fonte: analisemacro.com.br', size=8) +
  scale_fill_distiller(palette = "Greens", limits=c(0.60,0.85),
                       name="IDH")+
  theme_minimal()+
  no_axis


#Com Plotly

df <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/2011_us_ag_exports.csv")
df$hover <- with(df, paste(state, '<br>', "Beef", beef, "Dairy", dairy, "<br>",
                           "Fruits", total.fruits, "Veggies", total.veggies,
                           "<br>", "Wheat", wheat, "Corn", corn))
#Linha branca nas fronteiras
l <- list(color = toRGB("white"), width = 2)
# Opções
g <- list(
  scope = 'usa',
  projection = list(type = 'albers usa'),
  showlakes = TRUE,
  lakecolor = toRGB('white')
)

fig <- plot_geo(df, locationmode = 'USA-states')
fig <- fig %>% add_trace(
  z = ~total.exports, text = ~hover, locations = ~code,
  color = ~total.exports, colors = 'Purples'
)
fig <- fig %>% colorbar(title = "Millions USD")
fig <- fig %>% layout(
  title = '2011 US Agriculture Exports by State<br>(Hover for breakdown)',
  geo = g
)

fig

###Exemplo de mapa de bolha com Plotly

dfb<-world.cities[world.cities$country.etc=="Brazil",]
dfb$hover <- paste(dfb$name, "População", dfb$pop/1e6, "milhões")

dfb$q <- with(dfb, cut(pop, quantile(pop), include.lowest = T))
levels(dfb$q) <- paste(c("1st", "2nd", "3rd", "4th"), "Quantile")
dfb$q <- as.ordered(dfb$q)

g <- list(
  scope = 'south america',
  showland = TRUE,
  landcolor = toRGB("gray85"),
  subunitwidth = 1,
  countrywidth = 1,
  subunitcolor = toRGB("white"),
  countrycolor = toRGB("white")
)


fig <- plot_geo(dfb)
fig <- fig %>% add_trace(
  x = ~long, y = ~lat, text = ~hover,
  marker = list(size = sqrt(dfb$pop/10000) + 1, line = list(width = 0)),
  color = ~q, type = 'scattergeo'
)
fig <- fig %>% layout(title = 'Populações <BR> clicar nas legendas para selecionar', geo = ~g)

fig


#Outros assuntos: Shiny Apps, Dash e automação com Rscript

#Principais fOntes:
#https://r4ds.had.co.nz/data-visualisation.html
#https://cran.r-project.org/web/packages/geobr/vignettes/intro_to_geobr.html
#https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
#http://yulab-smu.top/treedata-book/chapter4.html
#https://plotly-r.com/overview.html
#https://rpubs.com/andreasancheztapia/176768
