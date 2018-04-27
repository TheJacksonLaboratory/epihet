# -------------------
# define our includes
# -------------------
OBJ_DIR = ../../../obj/
INCLUDES = -Isrc/ -Iinclude/

# ----------------------------------
# define our source and object files
# ----------------------------------

SOURCES_BAMS= src/api/BamAlignment.cpp \
         src/api/BamMultiReader.cpp \
         src/api/BamReader.cpp \
         src/api/BamWriter.cpp \
         src/api/SamHeader.cpp \
         src/api/SamProgram.cpp \
         src/api/SamProgramChain.cpp \
         src/api/SamReadGroup.cpp \
         src/api/SamReadGroupDictionary.cpp \
         src/api/SamSequence.cpp \
         src/api/SamSequenceDictionary.cpp \
         src/api/internal/bam/BamHeader_p.cpp \
         src/api/internal/bam/BamMultiReader_p.cpp \
         src/api/internal/bam/BamRandomAccessController_p.cpp \
         src/api/internal/bam/BamReader_p.cpp \
         src/api/internal/bam/BamWriter_p.cpp \
         src/api/internal/index/BamIndexFactory_p.cpp \
         src/api/internal/index/BamStandardIndex_p.cpp \
         src/api/internal/index/BamToolsIndex_p.cpp \
         src/api/internal/io/BamDeviceFactory_p.cpp \
         src/api/internal/io/BamFile_p.cpp \
         src/api/internal/io/BamFtp_p.cpp \
         src/api/internal/io/BamHttp_p.cpp \
         src/api/internal/io/BamPipe_p.cpp \
         src/api/internal/io/BgzfStream_p.cpp \
         src/api/internal/io/ByteArray_p.cpp \
         src/api/internal/io/HostAddress_p.cpp \
         src/api/internal/io/HostInfo_p.cpp \
         src/api/internal/io/HttpHeader_p.cpp \
         src/api/internal/io/ILocalIODevice_p.cpp \
         src/api/internal/io/RollingBuffer_p.cpp \
         src/api/internal/io/TcpSocketEngine_p.cpp \
         src/api/internal/io/TcpSocketEngine_unix_p.cpp \
         src/api/internal/io/TcpSocket_p.cpp \
         src/api/internal/sam/SamFormatParser_p.cpp \
         src/api/internal/sam/SamFormatPrinter_p.cpp \
         src/api/internal/sam/SamHeaderValidator_p.cpp \
         src/api/internal/utils/BamException_p.cpp

OBJECTS_BAMS= $(SOURCES_BAMS:.cpp=.o)
LIBRARY=lib/libbamtools.a
CXXFLAGS=-fPIC
all: $(LIBRARY)
    
.PHONY: all



$(LIBRARY): $(OBJECTS_BAMS)
	@echo " * linking $(LIBRARY)"
	ar cr $(LIBRARY) $(OBJECTS_BAMS)
	-ranlib $(LIBRARY)

src/api/BamAlignment.o: src/api/BamAlignment.cpp
	g++ -c src/api/BamAlignment.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/BamAlignment.o
src/api/BamMultiReader.o:
	g++ -c src/api/BamMultiReader.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/BamMultiReader.o
src/api/BamReader.o:
	g++ -c src/api/BamReader.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/BamReader.o
src/api/BamWriter.o:
	g++ -c src/api/BamWriter.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/BamWriter.o
src/api/SamProgram.o:
	g++ -c src/api/SamProgram.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamProgram.o
src/api/SamProgramChain.o:
	g++ -c src/api/SamProgramChain.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamProgramChain.o
src/api/SamReadGroup.o:
	g++ -c src/api/SamReadGroup.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamReadGroup.o
src/api/SamReadGroupDictionary.o:
	g++ -c src/api/SamReadGroupDictionary.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamReadGroupDictionary.o
src/api/SamSequence.o:
	g++ -c src/api/SamSequence.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamSequence.o
src/api/SamSequenceDictionary.o:
	g++ -c src/api/SamSequenceDictionary.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamSequenceDictionary.o
src/api/internal/bam/BamHeader_p.o:
	g++ -c src/api/internal/bam/BamHeader_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/bam/BamHeader_p.o
src/api/internal/bam/BamMultiReader_p.o:
	g++ -c src/api/internal/bam/BamMultiReader_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/bam/BamMultiReader_p.o
src/api/internal/bam/BamRandomAccessController_p.o:
	g++ -c src/api/internal/bam/BamRandomAccessController_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/bam/BamRandomAccessController_p.o
src/api/internal/bam/BamReader_p.o:
	g++ -c src/api/internal/bam/BamReader_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/bam/BamReader_p.o
src/api/internal/bam/BamWriter_p.o:
	g++ -c src/api/internal/bam/BamWriter_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/bam/BamWriter_p.o
src/api/internal/index/BamIndexFactory_p.o:
	g++ -c src/api/internal/index/BamIndexFactory_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/index/BamIndexFactory_p.o
src/api/internal/index/BamStandardIndex_p.o:
	g++ -c src/api/internal/index/BamStandardIndex_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/index/BamStandardIndex_p.o
src/api/internal/index/BamToolsIndex_p.o:
	g++ -c src/api/internal/index/BamToolsIndex_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/index/BamToolsIndex_p.o
src/api/internal/io/BamDeviceFactory_p.o:
	g++ -c src/api/internal/io/BamDeviceFactory_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/BamDeviceFactory_p.o
src/api/internal/io/BamFile_p.o:
	g++ -c src/api/internal/io/BamFile_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/BamFile_p.o
src/api/internal/io/BamFtp_p.o:
	g++ -c src/api/internal/io/BamFtp_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/BamFtp_p.o
src/api/internal/io/BamHttp_p.o:
	g++ -c src/api/internal/io/BamHttp_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/BamHttp_p.o
src/api/internal/io/BamPipe_p.o:
	g++ -c src/api/internal/io/BamPipe_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/BamPipe_p.o
src/api/internal/io/ByteArray_p.o:
	g++ -c src/api/internal/io/ByteArray_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/ByteArray_p.o
src/api/internal/io/HostAddress_p.o:
	g++ -c src/api/internal/io/HostAddress_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/HostAddress_p.o
src/api/internal/io/HostInfo_p.o:
	g++ -c src/api/internal/io/HostInfo_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/HostInfo_p.o
src/api/internal/io/HttpHeader_p.o:
	g++ -c src/api/internal/io/HttpHeader_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/HttpHeader_p.o
src/api/internal/io/ILocalIODevice_p.o:
	g++ -c src/api/internal/io/ILocalIODevice_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/ILocalIODevice_p.o
src/api/internal/io/RollingBuffer_p.o:
	g++ -c src/api/internal/io/RollingBuffer_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/RollingBuffer_p.o
src/api/internal/io/TcpSocketEngine_p.o:
	g++ -c src/api/internal/io/TcpSocketEngine_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/TcpSocketEngine_p.o
src/api/internal/io/TcpSocketEngine_unix_p.o:
	g++ -c src/api/internal/io/TcpSocketEngine_unix_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/TcpSocketEngine_unix_p.o
src/api/internal/io/TcpSocket_p.o:
	g++ -c src/api/internal/io/TcpSocket_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/TcpSocket_p.o
src/api/internal/sam/SamFormatParser_p.o:
	g++ -c src/api/internal/sam/SamFormatParser_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/sam/SamFormatParser_p.o
src/api/internal/sam/SamFormatPrinter_p.o:
	g++ -c src/api/internal/sam/SamFormatPrinter_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/sam/SamFormatPrinter_p.o
src/api/internal/sam/SamHeaderValidator_p.o:
	g++ -c src/api/internal/sam/SamHeaderValidator_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/sam/SamHeaderValidator_p.o
src/api/internal/utils/BamException_p.o:
	g++ -c src/api/internal/utils/BamException_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/utils/BamException_p.o
src/api/SamHeader.o:
	g++ -c src/api/SamHeader.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/SamHeader.o
src/api/internal/io/BgzfStream_p.o:
	g++ -c src/api/internal/io/BgzfStream_p.cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES) -o src/api/internal/io/BgzfStream_p.o


.PHONY: api
api:
	@echo " * Creating BamTools API"
	@mkdir -p lib
	@mkdir -p include
	@mkdir -p include/api
	@mkdir -p include/shared
	@cp src/api/*.h include/api
	@cp src/shared/*.h include/shared

.PHONY: clean_api
clean_api:
	@echo " * Cleaning-up"
	@rm -rf  ./src/api/*.o ./src/api/internal/*/*.o

