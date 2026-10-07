# Chaquopy: heavy reflection/JNI use, keep everything
-keep class com.chaquo.python.** { *; }

# App bridge classes (ZkBridge is the Java<->Python bridge)
-keep class com.witness.app.** { *; }
