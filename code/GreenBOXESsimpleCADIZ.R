# CAJAS VERDES ANDALUCIA
## Configuración y librerías

library(tidyverse)   # manipulación, análisis y visualización de datos (dplyr, ggplot2, etc.)
library(geosphere)   # cálculos geodésicos (distancias, áreas, ángulos, etc. sobre la superficie terrestre)
library(sf)          # Manejo y análisis de datos espaciales (shapefiles, geometrías, operaciones espaciales)
library(rstudioapi)  # Permite interactuar con RStudio (por ejemplo, obtener la ruta del script actual)
library(lubridate)   # Facilita el trabajo con fechas y horas (creación, conversión, extracción de componentes temporales)
library(tools)
library(ggplot2)   # visualización
library(rnaturalearth)       # mapas base
library(rnaturalearthdata)   # datos de Natural Earth
library(rnaturalearthhires) 
library(dplyr)

rm(list=ls(all=TRUE))# remove all previous functions

setwd(dirname(getActiveDocumentContext()$path))
getwd()

## FUNCIONES
create_cfr <- function(CODBU) {
  # Convertir a cadena
  CODBU <- as.character(CODBU)
  # Obtener la longitud del CODBU
  longitud <- nchar(CODBU)
  # Agregar ceros
  ceros <- paste0(rep("0", 9 - longitud), collapse = "")
  # Formatear el código completo
  # collapse="" une todos los elementos de un vector en un solo texto sin separarlos.
  # Con collapse="", convierte "0", "0", "0", "0" en "0000" como una sola cadena continua.
  codigo_formateado <- paste0("ESP", ceros, CODBU)
  # Asegurarse de que la longitud sea de 12 caracteres
  return(substr(codigo_formateado, 1, 12))
}






## Cargamos el fichero de cajas verdes
csv_entrada <- "jun22.csv"  #AQui ponemos el csv que vamos a procesar
data <- read_csv2(csv_entrada)
nomfichero <- file_path_sans_ext(csv_entrada)
# Leer el archivo CSV. Hay que comprobar que separadores de campos y de
# decimales tiene.
#
# Usaremos por convención el ; como separador de campos y el . como separador de
# decimales.
#

data <- unique (data) #Eliminamos duplicados
data <- data[, -c(4:6,11:14)] #BOrramos columnas no deseadas

#Renombrando columnas
data <- data %>% rename(SI_LATI = LATITUD, SI_LONG = LONGITUD, SI_SP = VELOCIDAD, SI_HE = RUMBO)
data <- data %>% rename(SI_DATE = FECHA, SI_TIME = HORA)


csv_formateado <- paste0(nomfichero,"formateado", ".csv")

write.table(data,
            file = csv_formateado,
            sep = ";",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# Trabajamos con el fichero formateado correctamente para evitar posteriores
# problemas.

vessel_track_df <- read_delim(csv_formateado, delim = ";")


# Vemos que la primera columna es el codigo de buque, hay que cambiarlo por
# el CFR. Aplicamos la función al dataframe

vessel_track_df$VE_REF <- sapply(vessel_track_df$CFPO, create_cfr)
vessel_track_df$VE_REF <- as.factor(vessel_track_df$VE_REF)

#Borro CFPO
vessel_track_df$CFPO<- NULL

# Creamos una columna de timestamp a partir de la fecha y la hora. Para ello
# usamos el paquete lubridate del tidyverse. OJO CON EL TIMEZONE.
#Un timestamp suele tener: Fecha: día, mes, año Hora: horas, minutos, segundos

vessel_track_df <- vessel_track_df %>%
  mutate(SI_TIMESTAMP = dmy_hms(paste(vessel_track_df$SI_DATE,
                                      vessel_track_df$SI_TIME),
                                tz = "Europe/Madrid")) 


#Ordenamos el dataset para visualizarlo
vessel_track_df <- vessel_track_df %>%
  arrange(VE_REF, SI_TIMESTAMP) %>%       # ordena las filas
  select(VE_REF, SI_TIMESTAMP, everything())  # mueve esas columnas al inicio

# Calculamos diferencias de tiempo y velocidades (en nudos)
# SI_TDIFF: diferencia de tiempo entre el timestamp actual y el anterior 
# SI_SPCA: calcula la velocidad = distancia / tiempo, Se multiplica por 1.94384 para convertir de m/s a nudos (1 m/s ≈ 1.94384 nudos).
# 1 m/s=3600/1852nudos=1.94384 nudos
# Un nudo (kt) se define como una milla náutica por hora.
# 1 milla náutica = 1852 m
# 1 hora = 3600 s
vessel_track_df <- vessel_track_df %>%
  group_by(VE_REF, floor_date(vessel_track_df$SI_TIMESTAMP, unit = "day")) %>%  #agrupa por barco y día
  mutate(
    SI_DISTANCECA = distHaversine(cbind(lag(SI_LONG), lag(SI_LATI)),  #Calcula la distancia en metros entre cada punto y el anterior, usando la fórmula de Haversine (para distancias geodésicas sobre la Tierra).
                                  cbind(SI_LONG, SI_LATI))
  ) %>%
  mutate(
    SI_TDIFF = (SI_TIMESTAMP - lag(SI_TIMESTAMP)),
    SI_SPCA = (SI_DISTANCECA / as.numeric(SI_TDIFF)) * 1.94384
  ) %>%
  ungroup()


#borramos la columna de floor_date
vessel_track_df <- vessel_track_df[, -9]


# Convertimos las columnas de latitud y longitud a geometrías POINT
# usando la librería sf.

points <- vessel_track_df %>%
  select(SI_LONG, SI_LATI) %>%
  st_as_sf(coords = c("SI_LONG", "SI_LATI"), crs = 4326) #El resultado es un objeto espacial coords con puntos geográficos (longitud-latitud).



# st_coordinates() extrae las coordenadas de los puntos y las devuelve como una matriz numérica.
coords_matrix <- st_coordinates(points) #Esta matriz se puede usar fácilmente para cálculos matemáticos o geométricos



# Calculamos el rumbo medio entre puntos consecutivos usando bearing
# ---------------------------------------------------------------
# La función bearing() nos da la dirección inicial para ir de un punto p1 a un punto p2
# sobre la superficie de la Tierra, siguiendo la ruta más corta (geodésico).
#
# Es decir, calcula el ángulo en grados entre los dos puntos considerando que la Tierra
# es curva (elipsoide WGS84), no plana.
#
# Aplicando esto a todos los puntos consecutivos de una trayectoria:
# - Se obtiene el rumbo del punto 1 al 2, del 2 al 3, del 3 al 4, y así sucesivamente.
# - Este rumbo indica hacia dónde "apunta" el movimiento al inicio de cada segmento.


bearings <- bearing(coords_matrix[-nrow(coords_matrix), ], coords_matrix[-1, ])


vessel_track_df$SI_COG <- c(NA, if_else(bearings < 0, bearings + 360, bearings))
#Course Over Ground (COG) “Rumbo sobre el suelo”.
#Si el rumbo (bearing) es negativo, le suma 360 para convertirlo a un ángulo positivo en el rango [0, 360).
#Si no es negativo, deja el valor tal cual.


# Eliminamos la primera fila ya que tendrá un NA creado al calcular la velocidad
vessel_track_df <- na.omit(vessel_track_df)

# Añadimos las columnas que es necesario completar y estimar. Estas por defecto
# estarán vacías El código FT_REF estará compuesto por 3 letras (abreviatura de
# centro) y números correlativos unidos por _
#
# EJEMPLO:
#
# SAN_0001 MUR_0004 BAL_0052 CAN_0456 etc...

# Marea
vessel_track_df$FT_REF <- NA
vessel_track_df$FT_REF <- as.factor(vessel_track_df$FT_REF)
#Creamos codigo de marea
vessel_track_df <- vessel_track_df %>%
  mutate(
    FT_REF = paste0(VE_REF, "_", format(as.Date(SI_TIMESTAMP, tz = "Europe/Madrid"), "%Y_%m_%d"))
  )



# Metier nivel 4
vessel_track_df$LE_MET4 <- NA
vessel_track_df$LE_MET4 <- as.factor(vessel_track_df$LE_MET4)

#Metier nivel 6
vessel_track_df$LE_MET6 <- NA
vessel_track_df$LE_MET6 <- as.factor(vessel_track_df$LE_MET6)

# Estado
#
# PESCA: TRUE
# NO PESCA: FALSE

vessel_track_df$SI_FSTATUS <- NA
vessel_track_df$SI_FSTATUS <- as.factor(vessel_track_df$SI_FSTATUS)

# Observador
#
# Observador embarcado : TRUE
# Observador no embarcado: FALSE

vessel_track_df$SU_ISOB <- NA
vessel_track_df$SU_ISOB <- as.factor(vessel_track_df$SU_ISOB)


# OGT
# Equipado: TRUE
# No equipado: FALSE
vessel_track_df$SI_OGT <- NA
vessel_track_df$SI_OGT <- as.factor(vessel_track_df$SI_OGT)


# En los GPS podríamos identificar también las operaciones de pesca
# En las cajas verdes este campo estará vacío, pero por homogeneidad
# de formato incluimos también esta columna
#
# LARGADA: SE
# VIRADA: HA
# ESPERA: WT (tiempo entre fin de largada e inicio de virada)
# NAVEGACION: ST
#
# Añadimos la culumna apropiada, que estará inicialmente vacía
vessel_track_df <- vessel_track_df %>%
  mutate(SI_FOPER = NA)

vessel_track_df$SI_FOPER <- as.factor(vessel_track_df$SI_FOPER)

# Nos quedamos con el dataframe como nos interesa

# vessel_track_df <- vessel_track_df %>%
#   select(VE_REF, FT_REF, SI_TIMESTAMP, SI_LATI, SI_LONG, SI_SP, SI_SPCA, SI_HE,
#          SI_COG, SI_DISTANCECA, SI_TDIFF, LE_MET4, LE_MET6,
#          SI_FSTATUS, SI_FOPER, SU_ISOB, SI_OGT)


vessel_track_df <- vessel_track_df %>%
select(VE_REF, FT_REF, SI_TIMESTAMP, SI_LATI, SI_LONG, SI_SP, SI_SPCA, SI_HE,
       SI_COG, SI_DISTANCECA, SI_TDIFF, LE_MET4, LE_MET6,
       SI_FSTATUS, SI_FOPER, SU_ISOB, SI_OGT) %>%
  mutate(
    # Convertimos timestamp a formato dd/mm/yyyy HH:MM:SS con / para Access
    SI_TIMESTAMP = format(SI_TIMESTAMP, "%d/%m/%Y %H:%M:%S")
  )


# Podemos guardar ya el fichero en formato csv
write.table(vessel_track_df,
            file = csv_formateado,
            sep = ";",
            dec = ",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)


# O en formato RDS (menos tamaño y 100% R, mantiene los tipos de las variables)
rds_formateado <- paste0(file_path_sans_ext(csv_formateado), ".rds")
saveRDS(vessel_track_df, file = rds_formateado)


# #Vamos a seleccionar puntos de cabo pino hacia algeciras
# names(vessel_track_df)
# barcos_filtrados <- vessel_track_df %>%  filter(SI_LONG < -4.74 & SI_LONG >-5.60)
# 
# write.table(barcos_filtrados,
#             file = "filtrados.csv",
#             sep = ";",
#             dec = ",",
#             row.names = FALSE,
#             col.names = TRUE,
#             quote = FALSE)





