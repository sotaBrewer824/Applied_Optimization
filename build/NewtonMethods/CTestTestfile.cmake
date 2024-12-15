# CMake generated Testfile for 
# Source directory: /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods
# Build directory: /Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/NewtonMethods
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(SpringElement2DWithLengthPSDHess.CheckHessian "/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/Build/bin/NewtonMethods-test" "--gtest_filter=SpringElement2DWithLengthPSDHess.CheckHessian" "--gtest_color=yes")
set_tests_properties(SpringElement2DWithLengthPSDHess.CheckHessian PROPERTIES  DEF_SOURCE_LINE "unit_tests.cc:124" SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.30.3/share/cmake/Modules/GoogleTest.cmake;438;add_test;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;28;gtest_add_tests;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;0;")
add_test(StandardNewton.CheckAlgorithm "/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/Build/bin/NewtonMethods-test" "--gtest_filter=StandardNewton.CheckAlgorithm" "--gtest_color=yes")
set_tests_properties(StandardNewton.CheckAlgorithm PROPERTIES  DEF_SOURCE_LINE "unit_tests.cc:144" SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.30.3/share/cmake/Modules/GoogleTest.cmake;438;add_test;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;28;gtest_add_tests;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;0;")
add_test(ProjectedNewton.CheckAlgorithmOnSpringElementWithoutLength "/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/Build/bin/NewtonMethods-test" "--gtest_filter=ProjectedNewton.CheckAlgorithmOnSpringElementWithoutLength" "--gtest_color=yes")
set_tests_properties(ProjectedNewton.CheckAlgorithmOnSpringElementWithoutLength PROPERTIES  DEF_SOURCE_LINE "unit_tests.cc:170" SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.30.3/share/cmake/Modules/GoogleTest.cmake;438;add_test;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;28;gtest_add_tests;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;0;")
add_test(ProjectedNewton.CheckAlgorithmOnSpringElementWithLength "/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/Build/bin/NewtonMethods-test" "--gtest_filter=ProjectedNewton.CheckAlgorithmOnSpringElementWithLength" "--gtest_color=yes")
set_tests_properties(ProjectedNewton.CheckAlgorithmOnSpringElementWithLength PROPERTIES  DEF_SOURCE_LINE "unit_tests.cc:193" SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.30.3/share/cmake/Modules/GoogleTest.cmake;438;add_test;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;28;gtest_add_tests;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;0;")
add_test(ProjectedNewton.CheckAlgorithmOnConstrainedSpringElement "/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/Build/bin/NewtonMethods-test" "--gtest_filter=ProjectedNewton.CheckAlgorithmOnConstrainedSpringElement" "--gtest_color=yes")
set_tests_properties(ProjectedNewton.CheckAlgorithmOnConstrainedSpringElement PROPERTIES  DEF_SOURCE_LINE "unit_tests.cc:215" SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.30.3/share/cmake/Modules/GoogleTest.cmake;438;add_test;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;28;gtest_add_tests;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;0;")
add_test(ProjectedNewton.CheckAlgorithmOnNonConvexFunction "/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/build/Build/bin/NewtonMethods-test" "--gtest_filter=ProjectedNewton.CheckAlgorithmOnNonConvexFunction" "--gtest_color=yes")
set_tests_properties(ProjectedNewton.CheckAlgorithmOnNonConvexFunction PROPERTIES  DEF_SOURCE_LINE "unit_tests.cc:239" SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]" _BACKTRACE_TRIPLES "/opt/homebrew/Cellar/cmake/3.30.3/share/cmake/Modules/GoogleTest.cmake;438;add_test;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;28;gtest_add_tests;/Users/xiazixuan/Desktop/CS_Self_Learning/Course_in_Bern/Applied_Optimizaton/Exercise11/aopt-exercise11/NewtonMethods/CMakeLists.txt;0;")
