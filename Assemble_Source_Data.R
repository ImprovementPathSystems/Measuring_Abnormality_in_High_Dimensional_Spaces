#' create a list of all the packages needed to run this R File. 
all_packages_needed <- c("reshape", "dplyr","proxy", "matrixStats","R.methodsS3", "zoo","gtools", "Hmisc","knitr","plyr")
#'
#'if an R package that is needed to run this code is not installed, install it. 
packages_to_install <- all_packages_needed[!(all_packages_needed %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) try(install.packages(packages_to_install),silent=TRUE)
#'
#'load all the packages needed. A FALSE value returned indicates the package failed to load. 
Load_Package_Silently <- function(Package_Name){suppressMessages(require(Package_Name,character.only = TRUE))}
sapply(all_packages_needed,Load_Package_Silently)#load all packages needed and don't print messages
#'
#'
#'Set variables needed to run this code
#Set the percent variance explained cutoff for TSQI, KmQI and KnQI. Value between 0 and 1.
# PCA_PercentVarExplainedCutoff <- .92

NormalSampleSize <- 32
Initial_Vector_Length <- 3500

NormalVectorFileDirectory <- read.csv('Source_File_Directory.csv', header=T, sep=',')
NormalVectorMatrix <- matrix(NA, nrow = Initial_Vector_Length, NormalSampleSize)

# ##Interlab File Directory
# InterlabDataFiles <- read.csv('Main_Source_Files/Interlab_FileDirectory_Final_Methodology.csv', header=T, sep=',')

# #load temporal spatial data
# refpop = as.data.frame(read.table('Main_Source_Files/TS_Normal_Variable_Matrix.csv', header=T, sep=','))
# refpopheights = as.data.frame(read.table('Main_Source_Files/TS_Normal_Height_Matrix_NoCadence_20141125.csv', header=T, sep=','))
# 



###############################################################################################################




#'Assemble all the normal vectors in a 3500 x NormalSampleSize matrix
for(n in 1:NormalSampleSize){
  
  CurrentNormVectorFilename <- NormalVectorFileDirectory$NormalVectorFileNames[n]
  CurrentNormData <- read.csv(paste('Source_Data/', CurrentNormVectorFilename,sep=''), header=T, sep=',')
  
  CurrentNormMean_AllVars <- subset(CurrentNormData, select=c(4))
  
  NormalVectorMatrix[,n] <- as.matrix(CurrentNormMean_AllVars)
}




#' #'**Begin the kinematic and kinetic calculations on the normal data**
#' ###############################################################################################################
#' #Calculate normal mean and std dev for all 3500 rows.
#' NormalMean <- as.matrix(rowMeans(NormalVectorMatrix))
#' NormalStdDev <- as.matrix(rowSds(NormalVectorMatrix))
#' 
#' 
#' ##Calculate z-scores on entire matrix
#' Z_Normal_AllVars <- apply(NormalVectorMatrix,2, function(x) (x-NormalMean)/NormalStdDev)






###############################################################################################################


#'###GDI Prep Calculations
##We need to calculate the previous methodololgies, the GDI and GPS to compare to the GQI
##This GDI section was created by Norman Dotson##

#######Setting Typically Developing Control Group Data###########
##################################
#Create list of measurements to be used
MeasurementList = c("L_Pelvis_Rotation"
                    ,"R_Pelvis_Fwd_Tilt"
                    ,"R_Pelvis_Lat_Tilt"
                    ,"L_HIP_Abd_ANG"
                    ,"L_HIP_Flex_ANG"
                    ,"L_HIP_Rot_ANG"
                    ,"L_KNEE_Flex_ANG"
                    ,"L_ANK_Flex_ANG"
                    ,"L_Foot_Orientation"
                    ,"R_HIP_Abd_ANG"
                    ,"R_HIP_Flex_ANG"
                    ,"R_HIP_Rot_ANG"
                    ,"R_KNEE_Flex_ANG"
                    ,"R_ANK_Flex_ANG"
                    ,"R_Foot_Orientation"
)

###########################
#######Variables###########
###########################
Vector_nrow = length(MeasurementList) * 100

#create an empty matrix with 100 * length(MeasurementAngleList)
#This will change according to vector variables selected
#Multiple Sample Size by 2 to account for left and right sides
GDI_GPS_NormalVectorMatrix <- matrix(NA, nrow = Vector_nrow, ncol=NormalSampleSize)


#Assemble all the normal vectors in a 3500 x NormalSampleSize matrix and assembler the GDI and GPS based on measurement list
for(n in 1:NormalSampleSize){
  
  CurrentNormVectorFilename <- NormalVectorFileDirectory$NormalVectorFileNames[n]
  CurrentNormData <- read.csv(paste('Source_Data/', CurrentNormVectorFilename,sep=''), header=T, sep=',')
  
  CurrentNormMean_AllVars <- subset(CurrentNormData, select=c(4))
  
  NormalVectorMatrix[,n] <- as.matrix(CurrentNormMean_AllVars)
  
  ###GDI and GPS Matrix##
  CurrentNormMean_ReqVars <- CurrentNormData[ which(CurrentNormData[,2] %in% MeasurementList),]
  ##Must include column 1 and 2 in 'order' to preserve a correct structure order
  #CurrentNormMean_ReqVars[,2] = as.character(CurrentNormMean_ReqVars[,2])
  CurrentNormMean_ReqVars <- CurrentNormMean_ReqVars[order(CurrentNormMean_ReqVars[,2],CurrentNormMean_ReqVars[,1]),]
  ##Pull the Mean values for each patient and append to NormalVectorMatrix
  GDI_GPS_NormalVectorMatrix[,n] <- as.matrix(CurrentNormMean_ReqVars[,4])
  
}

GDI_GPS_Var_Details <- CurrentNormMean_ReqVars[,1:3]

#' ##GDI Calculations for Normals###
#' NormMeans_ForGDIandGPS <- as.matrix(rowMeans(GDI_GPS_NormalVectorMatrix))
#' 
#' 
#' #'Able-Bodied Temporal Spatial Quality Index prep calculations
#' ###############################################################################################################
#' #Limit Temporal Spatial data based on number of normals in sample.
#' refpop <- refpop[which(refpop$NormID <= NormalSampleSize),]
#' refpopheights <- refpopheights[which(refpopheights$NormID <= NormalSampleSize),]
#' 
#' refpopheights <- subset(refpopheights, select=c(3,4,5,8,10,13,14))
#' 
#' refpop <- subset(refpop, select=c(3,4,5,8,10,13,14))
#' 
#' ##convert to a matrix
#' m_refpopheights <- data.matrix(refpopheights)
#' m_refpop <- data.matrix(refpop)
#' ##Normalize to height
#' refpopHeightNormalized <- m_refpopheights * m_refpop
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ###### Assemble Patient Data. ###### 
#' PatientVectorMatrix <- matrix(NA, nrow = Vector_nrow, ncol = length(InterlabDataFiles$VectorFileName))
#' PatientTemporalSpatialData <- matrix(NA, ncol= 7, nrow = length(InterlabDataFiles$VectorFileName))
#' 
#' #assemble interlab data into a single dataframe
#' #' ###Loop through all of the interlab data and calculate the GQI, GDI and GPS for each study.
#' #'These 10 subjects have been studied at 3 different time points over two days
#' #'This data is used to calculate the minimum detectable change (MDC) for the different methodologies
#' ###############################################################################################################
#' for(i in 1:length(InterlabDataFiles$VectorFileName)){
#'   
#'   #Read in interlab data
#'   CurrentVectorFilename <- InterlabDataFiles$VectorFileName[i]
#'   PatientData <- read.csv(paste('Interlab_Files/',CurrentVectorFilename,sep=''), header=T, sep=',')
#'   
#'   #Subset the data we need
#'   PatientMean_AllVars <- subset(PatientData, select=c(4))
#'   
#'   ###GDI and GPS Matrix##
#'   CurrentPatientMean_ReqVars <- PatientData[ which(PatientData[,2] %in% MeasurementList),]
#'   ##Must include column 1 and 2 in 'order' to preserve a correct structure order
#'   #CurrentNormMean_ReqVars[,2] = as.character(CurrentNormMean_ReqVars[,2])
#'   CurrentPatientMean_ReqVars <- CurrentPatientMean_ReqVars[order(CurrentPatientMean_ReqVars[,2],CurrentPatientMean_ReqVars[,1]),]
#'   PatientVectorMatrix[,i] <- CurrentPatientMean_ReqVars$Mean
#'   
#'   #Read in the patient's temporal Spatial data
#'   CurrentTemporalSpatialFilename <- InterlabDataFiles$TempSpatialFileName[i]
#'   CurrentTemporalSpatialData <- read.csv(paste('Interlab_Files/',CurrentTemporalSpatialFilename,sep=''), header=T, sep=',')
#'   
#'   #Read in the current height multiplier
#'   CurrentHeightMultiplier <- InterlabDataFiles$HeightMultiplier[i]
#'   
#'   CurrentHeightMultiplier <- as.numeric(as.character(CurrentHeightMultiplier))
#'   
#'   #Subset the data we need.
#'   CurrentTemporalSpatialData <- subset(CurrentTemporalSpatialData, select=c(1,2,3,6,8,11,12))
#'   
#'   ##Normalize to height
#'   CurrentTemporalSpatialData$L_Step_Length <- CurrentTemporalSpatialData$L_Step_Length * CurrentHeightMultiplier
#'   CurrentTemporalSpatialData$R_Step_Length <- CurrentTemporalSpatialData$R_Step_Length * CurrentHeightMultiplier
#'   CurrentTemporalSpatialData$L_Stride_Length <- CurrentTemporalSpatialData$L_Stride_Length * CurrentHeightMultiplier
#'   CurrentTemporalSpatialData$R_Stride_Length <- CurrentTemporalSpatialData$R_Stride_Length * CurrentHeightMultiplier
#'   CurrentTemporalSpatialData$L_Velocity <- CurrentTemporalSpatialData$L_Velocity * CurrentHeightMultiplier
#'   
#'   ##Per KK Edits, DO NOT NORMALIZE CADENCE TO HEIGHT.
#'   #CurrentTemporalSpatialData$L_Cadence <- CurrentTemporalSpatialData$L_Cadence * CurrentHeightMultiplier
#'   CurrentTemporalSpatialData$R_Step_Width <- CurrentTemporalSpatialData$R_Step_Width * CurrentHeightMultiplier
#'   
#'   PatientTemporalSpatialData[i,] <- as.matrix(CurrentTemporalSpatialData)
#' 
#' }