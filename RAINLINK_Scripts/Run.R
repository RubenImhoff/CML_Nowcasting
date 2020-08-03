## The RAINLINK package. Retrieval algorithm for rainfall mapping from microwave links 
## in a cellular communication network.
##
## Version 1.14
## Copyright (C) 2019 Aart Overeem
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## R program 'Run.R'.
## SCRIPT FOR RAINFALL ESTIMATION USING MICROWAVE LINKS.
## source("Run.R")
## Run this script by pasting the line above (source(...)) or by pasting parts of the script into the R shell.


# Note that it is not necessarily a problem if a function argument is not supplied to the function. If the
# function argument is not used, then there is no problem. Only be aware that you should use e.g.
# MaxFrequency=MaxFrequency. I.e. if you only supply MaxFrequency and the function argument before
# MaxFrequency is missing, than the function will not execute properly.


#############################################################
# 0. Load R libraries, parameter values, and other settings.#
# This also loads the RAINLINK package.                     #
#############################################################

source("Config.R") 

# If the time zone of the employed microwave link dataset is not the same as the (local) time zone used by R on your computer, set the time zone of the microwave link dataset:
# (this is important for functions RefLevelMinMaxRSL, WetDryNearbyLinkApMinMaxRSL and Interpolation):
#Sys.setenv(TZ='UTC')
# Otherwise RAINLINK can derive a wrong time interval length due to going to or from daylight saving time (DST). Timing of DST may be different between time zones, or one time zone may not have a change to/from DST.



############################
# 1. PreprocessingMinMaxRSL#
############################

# Load example data:
Linkdata <- read.table("../../4TUOpenData/4TU_upload/CMLs_20110609_20110911_21days.dat",header=TRUE)

# Add column with polarization if this column is not supplied in the link data:
if ("Polarization" %in% names(Linkdata)==FALSE)
{
   Linkdata$Polarization <- rep(NA,nrow(Linkdata))
}
# When no information on polarization is provided, the above code creates a column of NA for Polarization. In the function "RainRetrievalMinMaxRSL.R" links with
# NA values for polarization are processed with a & b values determined for vertically polarized signals.
# If information on polarization of links is available, use H for horizontally polarized & V for vertically polarized in Linkdata$Polarization.
# H, V & NA may occur in the same Linkdata file.

# Run R function:
StartTime <- proc.time()

DataPreprocessed <- PreprocessingMinMaxRSL(Data=Linkdata,MaxFrequency=MaxFrequency,MinFrequency=MinFrequency,verbose=TRUE)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))




############################################
# 2. WetDryNearbyLinkApMinMaxRSL (OPTIONAL)#
############################################

# Run R function:	
StartTime <- proc.time()

WetDry <- WetDryNearbyLinkApMinMaxRSL(Data=DataPreprocessed,CoorSystemInputData=NULL, 
MinHoursPmin=MinHoursPmin,PeriodHoursPmin=PeriodHoursPmin,Radius=Radius,Step8=Step8, 
ThresholdMedian=ThresholdMedian,ThresholdMedianL=ThresholdMedianL,ThresholdNumberLinks=ThresholdNumberLinks, 
ThresholdWetDry=ThresholdWetDry)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))




#######################
# 3. RefLevelMinMaxRSL#
#######################

# Run R function:
StartTime <- proc.time()

Pref <- RefLevelMinMaxRSL(Data=DataPreprocessed,Dry=WetDry$Dry,HoursRefLevel=HoursRefLevel,PeriodHoursRefLevel=PeriodHoursRefLevel)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))




#############################################################################################################
# 4. OutlierFilterMinMax (OPTIONAL) - Can only be applied when WetDryNearbyLinkApMinMaxRSL has been applied.#
#############################################################################################################

# Run R function:
DataOutlierFiltered <- OutlierFilterMinMaxRSL(Data=DataPreprocessed,F=WetDry$F,FilterThreshold=FilterThreshold)




######################
# 5. CorrectMinMaxRSL#
######################

# Run R function:
Pcor <- CorrectMinMaxRSL(Data=DataOutlierFiltered,Dry=WetDry$Dry,Pref=Pref)




############################
# 6. RainRetrievalMinMaxRSL#
############################

kRPowerLawDataH <- read.table(FileRainRetrHorizontal)
colnames(kRPowerLawDataH) <- c("f", "a", "b")

kRPowerLawDataV <- read.table(FileRainRetrVertical)
colnames(kRPowerLawDataV) <- c("f", "a", "b")


# Run R function:
StartTime <- proc.time()

Rmean <- RainRetrievalMinMaxRSL(Aa=Aa,alpha=alpha,Data=DataOutlierFiltered,kRPowerLawDataH=kRPowerLawDataH,kRPowerLawDataV=kRPowerLawDataV,PmaxCor=Pcor$PmaxCor,PminCor=Pcor$PminCor,Pref=Pref)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))



# Write path-average rainfall data to files:
# Duration of time interval of sampling strategy (min):
TIMESTEP <- 15	
	
# Location of output link data:
FolderRainEstimates <- paste("LinkPathRainDepths",TIMESTEP,"min",sep="")
ToFile = TRUE
if (ToFile)
{	
	# Create directory for output files:
	if(!dir.exists(FolderRainEstimates)){ dir.create(FolderRainEstimates) }
	# Write output to file
	ID <- unique(DataPreprocessed$ID)
	t <- sort(unique(DataPreprocessed$DateTime))
	t_sec <- as.numeric(as.POSIXct(as.character(t), format = "%Y%m%d%H%M"))
	dt <- min(diff(t_sec))
	
	for (i in 1 : length(t))
	{
		ind <- which(DataPreprocessed$DateTime == t[i])
		int_data <- data.frame(ID = DataPreprocessed$ID[ind], RainfallDepthPath = Rmean[ind] * dt / 3600, 
		PathLength = DataPreprocessed$PathLength[ind], XStart = DataPreprocessed$XStart[ind], 
		YStart = DataPreprocessed$YStart[ind], XEnd = DataPreprocessed$XEnd[ind], 
		YEnd = DataPreprocessed$YEnd[ind], IntervalNumber = rep(i, length(ind)), 
		Frequency = DataPreprocessed$Frequency[ind])
		
		Filename <- paste(FolderRainEstimates, "/linkdata_", t[i], ".dat", sep="")
		write.table(int_data, Filename, row.names = FALSE, col.names = TRUE, append = FALSE, quote = FALSE)
	}
}
# Note that the output files contain rainfall depths (mm). If these data are to be used for the interpolation, they must first be read ("Interpolation.R" does not read these files).
# Using the data for "Interpolation.R" requires a conversion from rainfall depth (mm) to rainfall intensity (mm/h).




###################
# 7. Interpolation#
###################

# Read grid onto which data are interpolated
RainGrid <- read.table(FileGrid, header = TRUE, sep=",")

# Duration of time interval of sampling strategy (min):
TIMESTEP <- 15		

# Location of output link data:
FolderRainMaps <- paste("RainMapsLinks",TIMESTEP,"min",sep="")

# Run R function:
StartTime <- proc.time()

Interpolation(Data=DataPreprocessed,CoorSystemInputData=NULL,idp=idp,IntpMethod=IntpMethod,nmax=nmax,
NUGGET=NUGGET,RANGE=RANGE,RainGrid=RainGrid,Rmean=Rmean,SILL=SILL,Variogram=Variogram,OutputDir=FolderRainMaps)

cat(sprintf("Finished. (%.1f seconds)\n",round((proc.time()-StartTime)[3],digits=1)))



