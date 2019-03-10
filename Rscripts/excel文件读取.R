两方面
1.是否安装JAVA，如果已经安装请检查JAVA是否符合R的版本。建议从新安装下JAVA：http://www.java.com/en/download/manual.jsp
2.不工作，在加载包之前，手动配置下java的位置
Sys.setenv(JAVA_HOME='C:\Program Files\Java\jre7') # for 64-bit version 
Sys.setenv(JAVA_HOME='C:\Program Files (x86)\Java\jre7') # for 32-bit version library(rJava)