#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=MinGW-Windows
CND_DLIB_EXT=dll
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Global_mod.o \
	${OBJECTDIR}/src/Inter_mod.o \
	${OBJECTDIR}/src/MATRIX.o \
	${OBJECTDIR}/src/Properties.o \
	${OBJECTDIR}/src/RWproc.o \
	${OBJECTDIR}/src/ReadWrite_f.o \
	${OBJECTDIR}/src/backsw.o \
	${OBJECTDIR}/src/bclq.o \
	${OBJECTDIR}/src/coeff.o \
	${OBJECTDIR}/src/cross.o \
	${OBJECTDIR}/src/dslucs_1.o \
	${OBJECTDIR}/src/output_old.o \
	${OBJECTDIR}/src/output_up.o \
	${OBJECTDIR}/src/rdbnds.o \
	${OBJECTDIR}/src/rdlito.o \
	${OBJECTDIR}/src/rivtox_dr.o \
	${OBJECTDIR}/src/rwBouss2DBinary.o \
	${OBJECTDIR}/src/rwIMMSPAscii.o \
	${OBJECTDIR}/src/rwIMMSPBinary.o \
	${OBJECTDIR}/src/rwSMSGenericBinary.o \
	${OBJECTDIR}/src/sedeq_cell.o \
	${OBJECTDIR}/src/solvl.o \
	${OBJECTDIR}/src/spch_bq.o \
	${OBJECTDIR}/src/srcsed.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/rivtoxclean.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/rivtoxclean.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/rivtoxclean ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/src/Global_mod.o: src/Global_mod.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/Global_mod.o src/Global_mod.f90

${OBJECTDIR}/src/Inter_mod.o: src/Inter_mod.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/Inter_mod.o src/Inter_mod.f90

${OBJECTDIR}/src/MATRIX.o: src/MATRIX.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MATRIX.o src/MATRIX.f

${OBJECTDIR}/src/Properties.o: src/Properties.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/Properties.o src/Properties.f90

${OBJECTDIR}/src/RWproc.o: src/RWproc.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/RWproc.o src/RWproc.f90

${OBJECTDIR}/src/ReadWrite_f.o: src/ReadWrite_f.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/ReadWrite_f.o src/ReadWrite_f.f90

${OBJECTDIR}/src/backsw.o: src/backsw.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/backsw.o src/backsw.f

${OBJECTDIR}/src/bclq.o: src/bclq.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/bclq.o src/bclq.f

${OBJECTDIR}/src/coeff.o: src/coeff.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/coeff.o src/coeff.f

${OBJECTDIR}/src/cross.o: src/cross.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/cross.o src/cross.f

${OBJECTDIR}/src/dslucs_1.o: src/dslucs_1.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/dslucs_1.o src/dslucs_1.f

${OBJECTDIR}/src/output_old.o: src/output_old.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/output_old.o src/output_old.f90

${OBJECTDIR}/src/output_up.o: src/output_up.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/output_up.o src/output_up.f90

${OBJECTDIR}/src/rdbnds.o: src/rdbnds.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rdbnds.o src/rdbnds.f90

${OBJECTDIR}/src/rdlito.o: src/rdlito.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rdlito.o src/rdlito.f

${OBJECTDIR}/src/rivtox_dr.o: src/rivtox_dr.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rivtox_dr.o src/rivtox_dr.f90

${OBJECTDIR}/src/rwBouss2DBinary.o: src/rwBouss2DBinary.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rwBouss2DBinary.o src/rwBouss2DBinary.f90

${OBJECTDIR}/src/rwIMMSPAscii.o: src/rwIMMSPAscii.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rwIMMSPAscii.o src/rwIMMSPAscii.f90

${OBJECTDIR}/src/rwIMMSPBinary.o: src/rwIMMSPBinary.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rwIMMSPBinary.o src/rwIMMSPBinary.f90

${OBJECTDIR}/src/rwSMSGenericBinary.o: src/rwSMSGenericBinary.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/rwSMSGenericBinary.o src/rwSMSGenericBinary.f90

${OBJECTDIR}/src/sedeq_cell.o: src/sedeq_cell.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/sedeq_cell.o src/sedeq_cell.f90

${OBJECTDIR}/src/solvl.o: src/solvl.f 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/solvl.o src/solvl.f

${OBJECTDIR}/src/spch_bq.o: src/spch_bq.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/spch_bq.o src/spch_bq.f90

${OBJECTDIR}/src/srcsed.o: src/srcsed.f90 
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/srcsed.o src/srcsed.f90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/rivtoxclean.exe
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
