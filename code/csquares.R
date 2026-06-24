library(sf)          # Manejo y análisis de datos espaciales (shapefiles, geometrías, operaciones espaciales)
library(lwgeom)# R version 4.0.2 para calculo de areas
library(rstudioapi)  # Permite interactuar con RStudio (por ejemplo, obtener la ruta del script actual)
library(ggplot2)   # visualización
library(rnaturalearth)       # mapas base
library(rnaturalearthdata)   # datos de Natural Earth
library(rnaturalearthhires) 



rm(list=ls(all=TRUE))# remove all previous functions
setwd(dirname(getActiveDocumentContext()$path))
getwd()


# Csquare Code Function  ################################################################################
# modificada de de la funcion Csquare de VMSTOOLS
csquare.cod<-function(lon, lat, degrees){
  if (length(lon) != length(lat)) 
    stop("length of longitude not equal to length of latitude")
  if (!degrees %in% c(10, 5, 1, 0.5, 0.1, 0.05, 0.01,0.005,0.001)) 
    stop("degrees specified not in range: c(10,5,1,0.5,0.1,0.05,0.01,0.005,0.001)")
  dims <- length(lon)
  quadrants <- array(NA, dim = c(5, 6, dims), dimnames = list(c("globalQuadrant","intmQuadrant1", "intmQuadrant2", "intmQuadrant3","intmQuadrant4"),c("quadrantDigit", "latDigit", "lonDigit","latRemain", "lonRemain", "code"),seq(1, dims, 1)))
  quadrants["globalQuadrant", "quadrantDigit",] <- 4 - (((2 * floor(1 + (lon/200))) - 1) * ((2 * floor(1 + (lat/200))) + 1))
  quadrants["globalQuadrant", "latDigit", ] <- floor(abs(lat)/10)
  quadrants["globalQuadrant", "lonDigit", ] <- floor(abs(lon)/10)
  quadrants["globalQuadrant", "latRemain", ] <- round(abs(lat) - (quadrants["globalQuadrant", "latDigit", ] * 10), 7)
  quadrants["globalQuadrant", "lonRemain", ] <- round(abs(lon) - (quadrants["globalQuadrant", "lonDigit", ] * 10), 7)
  quadrants["globalQuadrant", "code", ] <- quadrants["globalQuadrant", "quadrantDigit", ] * 1000 + quadrants["globalQuadrant", "latDigit", ] * 100 + quadrants["globalQuadrant", "lonDigit", ]
  quadrants["intmQuadrant1", "quadrantDigit", ] <- (2 * floor(quadrants["globalQuadrant", "latRemain", ] * 0.2)) + floor(quadrants["globalQuadrant", "lonRemain", ] * 0.2) + 1
  quadrants["intmQuadrant1", "latDigit", ] <- floor(quadrants["globalQuadrant", "latRemain", ])
  quadrants["intmQuadrant1", "lonDigit", ] <- floor(quadrants["globalQuadrant", "lonRemain", ])
  quadrants["intmQuadrant1", "latRemain", ] <- round((quadrants["globalQuadrant","latRemain", ] - quadrants["intmQuadrant1", "latDigit", ]) * 10, 7)
  quadrants["intmQuadrant1", "lonRemain", ] <- round((quadrants["globalQuadrant", "lonRemain", ] - quadrants["intmQuadrant1", "lonDigit", ]) * 10, 7)
  quadrants["intmQuadrant1", "code", ] <- quadrants["intmQuadrant1", "quadrantDigit", ] * 100 + quadrants["intmQuadrant1", "latDigit", ] * 10 + quadrants["intmQuadrant1", "lonDigit", ]
  quadrants["intmQuadrant2", "quadrantDigit", ] <- (2 * floor(quadrants["intmQuadrant1", "latRemain", ] * 0.2)) + floor(quadrants["intmQuadrant1", "lonRemain", ] * 0.2) + 1
  quadrants["intmQuadrant2", "latDigit", ] <- floor(quadrants["intmQuadrant1", "latRemain", ])
  quadrants["intmQuadrant2", "lonDigit", ] <- floor(quadrants["intmQuadrant1", "lonRemain", ])
  quadrants["intmQuadrant2", "latRemain", ] <- round((quadrants["intmQuadrant1", "latRemain", ] - quadrants["intmQuadrant2", "latDigit", ]) * 10, 7)
  quadrants["intmQuadrant2", "lonRemain", ] <- round((quadrants["intmQuadrant1", "lonRemain", ] - quadrants["intmQuadrant2", "lonDigit", ]) * 10, 7)
  quadrants["intmQuadrant2", "code", ] <- quadrants["intmQuadrant2", "quadrantDigit", ] * 100 + quadrants["intmQuadrant2", "latDigit", ] * 10 + quadrants["intmQuadrant2", "lonDigit", ]
  quadrants["intmQuadrant3", "quadrantDigit", ] <- (2 * floor(quadrants["intmQuadrant2", "latRemain", ] * 0.2)) + floor(quadrants["intmQuadrant2", "lonRemain", ] * 0.2) + 1
  quadrants["intmQuadrant3", "latDigit", ] <- floor(quadrants["intmQuadrant2", "latRemain", ])
  quadrants["intmQuadrant3", "lonDigit", ] <- floor(quadrants["intmQuadrant2", "lonRemain", ])
  quadrants["intmQuadrant3", "latRemain", ] <- round((quadrants["intmQuadrant2","latRemain", ] - quadrants["intmQuadrant3", "latDigit", ]) * 10, 7)
  quadrants["intmQuadrant3", "lonRemain", ] <- round((quadrants["intmQuadrant2", "lonRemain", ] - quadrants["intmQuadrant3", "lonDigit", ]) * 10, 7)
  quadrants["intmQuadrant3", "code", ] <- quadrants["intmQuadrant3", "quadrantDigit", ] * 100 + quadrants["intmQuadrant3", "latDigit", ] * 10 + quadrants["intmQuadrant3", "lonDigit", ]
  quadrants["intmQuadrant4", "quadrantDigit", ] <- (2 * floor(quadrants["intmQuadrant3", "latRemain", ] * 0.2)) + floor(quadrants["intmQuadrant3", "lonRemain", ] * 0.2) + 1
  quadrants["intmQuadrant4", "latDigit", ] <- floor(quadrants["intmQuadrant3", "latRemain", ])
  quadrants["intmQuadrant4", "lonDigit", ] <- floor(quadrants["intmQuadrant3", "lonRemain", ])
  quadrants["intmQuadrant4", "latRemain", ] <- round((quadrants["intmQuadrant3","latRemain", ] - quadrants["intmQuadrant4", "latDigit", ]) * 10, 7)
  quadrants["intmQuadrant4", "lonRemain", ] <- round((quadrants["intmQuadrant3", "lonRemain", ] - quadrants["intmQuadrant4", "lonDigit", ]) * 10, 7)
  quadrants["intmQuadrant4", "code", ] <- quadrants["intmQuadrant4","quadrantDigit", ] * 100 + quadrants["intmQuadrant4","latDigit", ] * 10 + quadrants["intmQuadrant4", "lonDigit", ]
  if (degrees == 10) 
    CSquareCodes <- quadrants["globalQuadrant", "code",]
  if (degrees == 5) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "quadrantDigit", ], sep = "")
  if (degrees == 1) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], sep = "")
  if (degrees == 0.5) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], ":", quadrants["intmQuadrant2", "quadrantDigit", ], sep = "")
  if (degrees == 0.1) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], ":", quadrants["intmQuadrant2", "code", ], sep = "")
  if (degrees == 0.05) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], ":", quadrants["intmQuadrant2", "code", ], ":", quadrants["intmQuadrant3", "quadrantDigit", ], sep = "")
  if (degrees == 0.01) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], ":", quadrants["intmQuadrant2", "code", ], ":", quadrants["intmQuadrant3", "code", ], sep = "")
  if (degrees == 0.005) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], ":", quadrants["intmQuadrant2", "code", ], ":", quadrants["intmQuadrant3", "code", ], ":", quadrants["intmQuadrant4", "quadrantDigit", ], sep = "")
  if (degrees == 0.001) 
    CSquareCodes <- paste(quadrants["globalQuadrant", "code", ], ":", quadrants["intmQuadrant1", "code", ], ":", quadrants["intmQuadrant2", "code", ], ":", quadrants["intmQuadrant3", "code", ], ":", quadrants["intmQuadrant4", "code", ], sep = "")
  
  return(CSquareCodes)
}


# CFija directorio trabajo ################################################################################



resolution<-0.05 # millas nauticas de lado=resolution*60
# Para que este anidado y estandarizado 10,5,1,0.5,0.1,0.05,0.01,0.005,0.001)

#1 grado de latitud ≈ 60 millas náuticas, por definición (1 milla náutica = 1 minuto)
# Conversión de grados a millas náuticas y kilómetros:
# 1 grado = 60 millas náuticas (por definición)
# 0.05° × 60 = 3 millas náuticas
# 3 × 1.852 = 5.556 km
# → Cada lado del cuadrado del grid mide aproximadamente 5.56 km

print(paste0("Resolucion: Milas Nauticas: ",resolution*60,"  // Kilometros: ", resolution*60*1.852))
      

      
      

#---- Esquinas de la zona de estudio
#-esquina inferior izq 
x.ini.izq.inf<-(-5)#usar numeros enteros
y.ini.izq.inf<-(35)
#-esquina sup derecha
x.fin.dch.sup<-(-7)
y.fin.dch.sup<-(37)

xncell<-abs(x.ini.izq.inf-x.fin.dch.sup)/resolution
yncell<-abs(y.ini.izq.inf-y.fin.dch.sup)/resolution

df<-data.frame(ID=1:(xncell*yncell))

gr<-st_sfc(st_make_grid(cellsize=resolution,
                        offset=c(x.ini.izq.inf,y.ini.izq.inf),n=c(xncell,yncell),
                        what="polygons",square=T,crs=4326))

st_geometry(df)<-gr

#---- Add centroides
df$centroid <- st_centroid(df$geometry)
head(df)
df$LONG_CENT<-as.numeric(st_coordinates(st_centroid(df$geometry))[,1])
df$LAT_CENT<-as.numeric(st_coordinates(st_centroid(df$geometry))[,2])

#---- Add area de cada tesela
df$area_km2 <- as.numeric(st_area(df))/1000000# en km2

head(df)

#---- asigna Csquare
df$csq<-as.character(csquare.cod(df$LONG_CENT,df$LAT_CENT,resolution))

head(df)
#--- Guarda grid con datos en shapefiles
st_write(df, "migridcadiz.shp",delete_dsn = TRUE)



# Mapa de España
spain <- ne_countries(country = "Spain", scale = "medium", returnclass = "sf")

# Limites aproximados del sur de España
xlim <- c(-6, 0)    # longitudes: desde Huelva/Cádiz hasta Málaga/Almería
ylim <- c(36, 38.5) # latitudes: desde Cádiz hasta Granada/Almería

ggplot() +
  geom_sf(data = spain, fill = "gray95", color = "black", size = 0.2) + 
  geom_sf(data = df, fill = NA, color = "darkblue", size = 0.1) +
  coord_sf(xlim = xlim, ylim = ylim) +   # zoom a la zona sur
  labs(title = "Grid C-square (0.05°) - Sur de España") +
  theme_minimal()
