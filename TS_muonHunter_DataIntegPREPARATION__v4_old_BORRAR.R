# Datos recogidos con el MuonHunter
# This script is to read a csv file containing the data directly exported by the detector
#   using the python script export.py and prepare the data for
#   its analysis. The analysis is performed by other script
#   The Analysis is check for periodicity and generate a
#   time serie (with cph) (muon "Channel" (coincidences)) and a stl.
#
# This script also detects malfunctions and errors in data set
#
# INSTRUCTIONS
# Change the csv file name
# Change the values of lowlimit and highlimit. Values to filter the wrong cph
lowlimit <- 0
highlimit <- 1000
# Esta version incluye missing_s field in aux_count data frame
# It will use aggregate function to integrate seconds into minutes and hours
# ESTA A MEDIAS

# DESCRIPTION OF MAIN INPUT, OUTPUT & AUX DATA STRUCTURES USED
# (IN) muonDataFile: stores the info from the detector csv file
# (OUT) muon: data frame ($time, $cpm, $cph) 
#		minute by minute from the beginning to the end
#		stored in the file muonhunter-data-perMinute.csv
# (OUT) muon_h: data frame ($time $cph) generated from muon,
#		integrating the 60 mins. in a hour
#		stored in the file muonhunter-data-perHour.csv
# aux_count: auxiliar data frame. Used to calculate the 
#		counts in a second from the Total Hits column of
#		the detector csv file. To use it in integration
# 		in minutes (muon) and hours (muon_h)
#		Is a data frame (time, totalCount, periodCount)
#		$time: POSIXct, second recorded in detector's csv file
#		$totalCount: Total Hits column of detector's csv file
#			Since detector was started. Warning! reseted if 
#			detector is reseted. And some times wrong value
#		$periodCount: calculated from totalCount as increment
#			since previour record
# NEW
#		$missing_s: number of missing seconds since the previous
#			record. Usually 0, but some times 1, less times some
#			seconds, and in special situations 
#			high (db off or detector off)

# The script autodetect nrow = (records-1)/3 -1 by the header. /3 because of GM1, GM2, Muon
# And it changes start and end

# Load data
setwd("C:/Users/h.rio/OneDrive - GSD Educación/muonhunter_R")
# Modified on 2018/07/20 Always the same name for input data: muonhunter-data.csv
# muonDataFile<-read.table("muonhunter-data-prueba2.csv",
muonDataFile<-read.table("muonhunter-data.csv",
                  header=T, dec=".", sep=",",stringsAsFactors=F)
cat("File Lines Number", dim(muonDataFile)[1], "\n")

# How to calculate number of registers with Channel=Muon. 
#In a previous version this is used to generate a matrix
#registros <- dim(muonh)[1]/3
#registros

# muonh$Log.time[1] is the time for the first entry
cat("Starting time: ", muonDataFile$Log.time[1], "\n" )
# muonh$Log.time[dim(muonh)[1]] is the time for the last file entry
cat("End time: ", muonDataFile$Log.time[dim(muonDataFile)[1]], "\n")


# Data frame diferente del script de DataPreview
# Se eliminan las columnas gm1 y gm2
# para reducir la cantidad de memoria
# time is POSIXct date and time by minutes
# cpm counts in this minute
# cph counts acumulated in the last 60 minutes (last hour),
#   included the currently processed minute
# A MEDIAS
# 	missing_s number of missing seconds in the minute (discard if bad)
first_minute = as.POSIXct(muonDataFile$Log.time[1],format="%Y-%m-%d %H:%M", tz="GMT")
muon <- data.frame(
  time=seq(from = first_minute,
           to = as.POSIXct(muonDataFile$Log.time[dim(muonDataFile)[1]],format="%Y-%m-%d %H:%M", tz="GMT"), 
           by = "min"),
  cpm=0,
  cph=0,
  missing_s = 0)

###########################################################
# Code to process to calculate aux_count:
#		Calculate counts per period(second) from total hits
#		Identify missing seconds and other lacks in data
###########################################################
  
# Initialization of aux_count and auxiliar variables
# Fill the entries (by time) which have Muon (coincidence) data
aux_count = data.frame(
  time=muonDataFile[which(muonDataFile$Channel=="Muon"),3],
  totalCount=muonDataFile[which(muonDataFile$Channel=="Muon"),6],
  periodCount=NA,
  missing_s = -1
)
wrongCount_threshold = 20
missing_s_threshold = 10
missing_s_bw_entries_threshold = 30
prev_count = aux_count$totalCount[1]
aux_count$time <- as.POSIXlt(aux_count$time, tz = "GMT")
aux_count$periodCount[1] = 0
aux_count$missing_s[1] = 0


for (i in  2:length(aux_count$periodCount)) {
  missing_s_bw_entries = difftime(aux_count$time[i]-1, 
                                    aux_count$time[i-1], 
                                    units = "secs")
	if(missing_s_bw_entries > missing_s_bw_entries_threshold) {
	cat("More than", missing_s_bw_entries_threshold," s missing at ")
	print(aux_count$time[i])
	cat("since ")
	print(aux_count$time[i-1])
  }
  aux_count$missing_s[i] = missing_s_bw_entries
 # missing seconds stored in i-th entry corresponds to missing 
 #		seconds BEFORE the currently processed
  
  if (aux_count$missing_s[i] < missing_s_threshold)
  {		# There aren´t missing seconds
		#	or there are few seconds missing
	if ( abs(aux_count$totalCount[i] - prev_count) <= wrongCount_threshold )
	{
# Normal Behaviour
	    aux_count$periodCount[i] = aux_count$totalCount[i] - prev_count
		prev_count = aux_count$totalCount[i]
	}
	else{
# totalCount either too high or too small
# abs(aux_count$totalCount[i] - prev_count) > wrongCount_threshold
		cat("i ",i," totalCount ",aux_count$totalCount[i])
		cat(" prev_count ", prev_count)
		cat(" missing sec=", aux_count$missing_s[i])
		print(aux_count$time[i])
# Usually aux_count$totalCount[i+1] == aux_count$totalCount[i-1]
# 	But condition is modified to allow few counts
		if ( (aux_count$totalCount[i+1] - aux_count$totalCount[i-1]) < 4 )
		{
# Wrong detector data this second. prev_count is bad value.
#Remove this periodCount and fix the previous one
			aux_count$periodCount[i] = 0
			# prev_count is the i-1 value      
			prev_count = aux_count$totalCount[i-1]
			# wrong totalCount is left as error remainder
		}
		else {
# A big difference in counts should happen Only if there is not data recording.
# Then there are missing seconds and this part of the code is not reached
			cat("i=",i)
			cat(" aux_count$totalCount[", i+1, "]=", aux_count$totalCount[i+1])
			cat(" aux_count$totalCount[", i-1, "]=", aux_count$totalCount[i-1])
			cat(" Values differ in more than 4 counts")
			cat("PORQUE HA LLEGADO AQUI? NO DEBERIA\n")
		}
	}
 } # If there is not missing seconds
 else {
# THIS CODE is executed when no data are recorded in a while and the detector
#	continues working (on). 
# 	ejemplo datos de 16 de abril de 07:12:05 pasa a 07:53:07
#   Que hacer con las cuentas??? Repartir uniformemente en los segundos pasados???
#     1 segundo cada minuto con cuentas acumuladas?????

    # reset prev_count variable for muonhunter either reset or restart
    # Then start in 0 or low value and it maintains.
    # There is no incremental count and restart the prev_count variable
    aux_count$periodCount[i] = 0
    prev_count = 0
 }
} # for

# For testing
write.table(aux_count, file="aux_count.csv",sep = ",", row.names = FALSE)



############ MODIFICAR USANDO aggregate


i = 1 # index of aux_count (seconds)
# j index of muon vector (minutes)


# suma_cpm sumatorio de los primeros minutos para hacer la interpolación
#   de cph en la primera hora
suma_cpm = 0

# initialization of indexes when there are missing data (i.e. detector stopped)
j_last_rec = 0  # last minute before the missing period
j_new_rec = 0   # fisrt minute after the missing period

#
# First minute. Usually few seconds with data. Then, not usefull
size_aux_count <- length(aux_count$time)
while ( (i <= size_aux_count) &&
        ( as.POSIXct(aux_count$time[i]) < (first_minute+60) ) )
  {
    muon$cpm[1] = muon$cpm[1] + aux_count$periodCount[i]
	muon$missing_s[1] = muon$missing_s[1] + aux_count$missing_s[i]
	i = i + 1
  }
  #Interpolation of cph from the 1st minute
  muon$cph[1] = round(60*muon$cpm[1])

# j=2 is second minute
j = 2
# loop for minutes
while (j <= length(muon$time)){

#To stop in debugging
#if (j == 218)  {
#  cat("To stop\n")
#}
   # If there is one or more minutes missing BEFORE the processed one
   # Should be checked before the normal behaviour 

  if (aux_count$missing_s[i] > 60)
  {
  # At this point aux_count$time[i] inside the first minute of resuming data recording
  # but should be muon$time[j-1]+60
  # the entry should be stored in j + trunc(muon$missing_s[j]/60) is partially recorded (database registering stopped)
	# j_last_rec is the j index of the last minute with some records
	j_last_rec = j - 1
	# j_new_rec is the j index of the minute where records are newly resumed. Partial data
	j_new_rec = j_last_rec + trunc(aux_count$missing_s[i]/60)
  # Since j - trunc(muon$missing_s[j]/60)) +1 to j-1 are missing
  # j is partially recorded (database register restarted)
  #
   # set values for completely missing minutes
	for(k in ((j_last_rec + 1):(j_new_rec -1) ) ){
	  muon$cpm[k] = NA
	  muon$missing_s[k] = 60
	  muon$cph[k] = NA
	}
  # values for j_last_rec (last minute partially recorded). 
  # muon$missing_s[j_last_rec] = lo que ponga de perdidos antes + los que van hasta el final del minuto
  # IS TOO COMPLEX   muon$missing_s[j_last_rec] = muon$missing_s[j_last_rec] + sec_to_min_end
	muon$missing_s[j_last_rec] = 66
	muon$cph[j_last_rec] = muon$cph[j_last_rec-1] # simplificacion
	# muon$cpm calculado en el while
	
  # set correct values in j_new_rec entry IS TOO COMPLEX
  # muon$missing_s[j] = muon$missing_s[j] - 60 * trunc(muon$missing_s[j]/60) + 1
  # todavía hay que quitar los del último minuto antes de parar
	muon$missing_s[j_new_rec] = 66
  # To avoid NA when accumulating counts in the j_new_rec minute
	muon$cpm[j_new_rec] = 0
	# To give a value to cpm to avoid NA values after minutes without values
	# use last calculated value if passed less than an hour since last minute with records.
	# else, restart to 0
#	if  ((muon$time[j_last_rec]+60 > muon$time[j_new_rec]) || is.na(muon$cph[j_last_rec-1]) ) {
#		muon$cph[j_new_rec] = 0
#	}
#	else {
#		muon$cph[j_last_rec] = muon$cph[j_last_rec-1]
#	}

	# como muon$cph[j_new_rec-1] es NA cuando haga el if a continuación se queda NA
	# y se deja j para ese valor para que calcule j_new_rec en la siguiente iteracion (tendrá missing_s mayor de 66)	
	j = j_new_rec
	
	# entrada i con missing_s > 60 already processed. Pass to the next one
	i = i+1
  } # if aux_count$missing_s[i] > 60
 # loop to go across the aux_count vector (seconds)
 # grouping seconds in minutes
 while ( (i <= length(aux_count$time)) &&
	( as.POSIXct(aux_count$time[i]) < (muon$time[j]+60) ) )
  {
  
#To stop at a fixed time in debugging
aux_t <- strptime("2018-04-16 07:12:05", format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
if ( difftime(aux_count$time[i], aux_t, units = "secs") == 0) {
        cat("To stop\n")
}
  # NORMAL BEHAVIOUR
  # loop for seconds inside the minute
	muon$cpm[j] = muon$cpm[j] + aux_count$periodCount[i]
	muon$missing_s[j] = muon$missing_s[j] + aux_count$missing_s[i]
  
  i = i+1
  
 } # while inside aux_count data set
 	# cph adding the LAST 60 minutes
	#   should be the same as the detector registrates in minutes???
	#		creo que no por la diferencia entre tiempo de detector y tiempo S.O. raspberry
	if (j > 60)
	{ # One hour of data recording
		if (j_new_rec == j) { 
		# calculating cph after missing minutes
		# Minutes without data have cpm=NA (it is not the same than cpm=0 a minute without detections)
		# Assumption cpm=NA managed as having cpm=0 REVISAR ESTA HIPOTESIS p.ej. restar la última detección muon$cpm[j-61]
			if ( is.na(muon$cpm[j-60]) ) { # more than hour without data. cph is restarted 
				muon$cph[j] = muon$cpm[j]
			}	
			else { # less than 1 hour missing
				# start with the last cph recorded, remove minute 60 minutes ago and add new one
				# 	in the case of more than 1 hour missing, cpm=NA and was maneged in previous if
				muon$cph[j]=muon$cph[j_last_rec]-muon$cpm[j-60]+muon$cpm[j]
			}
		}
		else {
			if ( is.na(muon$cpm[j-60]) ) { # cpm = NA, trying to remove a non existing data. Assume it aa 0 
				muon$cph[j] = muon$cph[j-1] + muon$cpm[j]
			}
			else {
				# NORMALLY: remove the old minute (60 seconds ago) and add the new one
				muon$cph[j]=muon$cph[j-1]-muon$cpm[j-60]+muon$cpm[j]
			}
		}
	} # if 1st hour already passed
	else { #Before the first hour: Interpolate
		suma_cpm = suma_cpm + muon$cpm[j]
		muon$cph[j] = suma_cpm + round((60-j)*suma_cpm/j)
	} # else
 j = j + 1
} # while j (minutes)

#Saving 
cat("Saving file integrating from seconds to minutes\n")
write.table(muon,file="muonhunter-data-perMinute.csv", dec=".",
            sep=",", row.names = FALSE)


hist(muon$cpm)
hist(muon$missing_s)
hist(muon$cph)

first_hour = as.POSIXct(muonDataFile$Log.time[1],format="%Y-%m-%d %H", tz="GMT")
muon_h <- data.frame(
  time=seq(from = first_hour,
           to = as.POSIXct(muonDataFile$Log.time[dim(muonDataFile)[1]],format="%Y-%m-%d %H", tz="GMT"), 
           by = "hour"),
  cph=0,
  missing_s = 0)
  

i = 1 # index of muon data frame vectors (minutes)
# j index of muon_h data frame vectors (hours)
#
# Process the first hour. Usually without the 60 minutes
size_muon <-length(muon$time)
while ( (i <= size_muon) &&
        ( as.POSIXct(muon$time[i]) < (first_hour+3600) ) )
  {
    muon_h$cph[1] = muon_h$cph[1] + muon$cpm[i]
	muon_h$missing_s[1] = muon_h$missing_s[1] + muon$missing_s[i]
	i = i + 1
}
# Add the seconds before data recording

for (j in 2:length(muon_h$time)){
  # Arreglar bucle sin fin
  while ( (i <= size_muon) &&
          ( as.POSIXct(muon$time[i]) < (muon_h$time[j]+3600) ) )
  {
	if (is.na(muon$cpm[i])) 	{
		# if cpm=NA then is not added. Data will be not analysed if missing_s is too high
		muon_h$cph[j] = muon_h$cph[j] + 0
	}
	else {
		muon_h$cph[j] = muon_h$cph[j] + muon$cpm[i]
	}
	muon_h$missing_s[j] = muon_h$missing_s[j] + muon$missing_s[i]
    i = i + 1
  }

}

hist(muon_h$cph)

#Saving 
cat("Saving file integrating from seconds to hour\n")
write.table(muon_h,file="muonhunter-data-perHour.csv", dec=".",
            sep=",", row.names = FALSE)
#####################################################
# Acabar aqui preparacion de datos
#####################################################

