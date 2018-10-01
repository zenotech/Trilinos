INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING
                    ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GeneralizedAlpha
  SOURCES Tempus_GeneralizedAlphaTest.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
  TESTONLYLIBS tempus_test_models
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(Test_GeneralizedAlpha_CopyFiles
  DEST_FILES Tempus_GeneralizedAlpha_HarmonicOscillator_Damped_SecondOrder.xml
  EXEDEPS GeneralizedAlpha
  )
