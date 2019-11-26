cmake_minimum_required(VERSION 3.10)
project(ShapeProducer)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g")

# Cint Reflex genmap Cintex Rint Thread New rootcint rlibmap Core rmkdepend MathCore MathMore Matrix Minuit Minuit2 Fumili Physics MLP Quadp Foam Smatrix SPlot GenVector Genetic Unuran Hist HistPainter Spectrum SpectrumPainter Hbook Tree TreePlayer TreeViewer RIO SQLIO XMLIO XMLParser Net RootAuth Krb5Auth Netx SrvAuth RDAVIX NetxNG RHTTP Gpad Graf Postscript mathtext GX11 GX11TTF ASImage ASImageGui FITSIO Graf3d X3d Eve RGL GLEW FTGL Gviz3d Gui Ged FitPanel GuiBld GuiHtml Recorder SessionViewer Proof ProofPlayer ProofDraw ProofBench proofd xproofd XrdProofd Proofx Html EG VMC Geom GeomBuilder GeomPainter Gdml root minicern MemStat Table rootn.exe roots.exe ssh2rpd root.exe proofserv.exe hadd g2root h2root PyROOT TMVA RooFitCore RooFit RooStats HistFactory hist2workspace
find_package(ROOT REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(jitted_code_dir ${CMAKE_BINARY_DIR}/jitted_code)
add_subdirectory(jsoncpp)

# Set the include paths for the compiler
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(ShapeProducer SHARED src/InputManager.cxx src/Cutstring.cxx)
target_include_directories(ShapeProducer PRIVATE inc)
target_link_libraries(ShapeProducer ${ROOT_LIBRARIES} jsoncpp)

#add_executable(gamma gamma.C)

file(COPY shape_producer DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# Copy C++ files with extension .h and .cxx in the same directory
# in order not to use weird include paths when we use jit code
# in Python
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/inc/ ${jitted_code_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/src/ ${jitted_code_dir})
