organization := "com.github.biopet"
organizationName := "Biowdl"
name := "tasks"

biopetUrlName := "tasks"

startYear := Some(2018)

biopetIsTool := false
biopetIsPipeline := true

developers += Developer(id = "ffinfo",
                        name = "Peter van 't Hof",
                        email = "pjrvanthof@gmail.com",
                        url = url("https://github.com/ffinfo"))
developers += Developer(id = "rhpvorderman",
                        name = "Ruben Vorderman",
                        email = "r.h.p.vorderman@lumc.nl",
                        url = url("https://github.com/rhpvorderman"))

scalaVersion := "2.11.12"

libraryDependencies += "com.github.biopet" %% "biowdl-test-utils" % "0.2-SNAPSHOT" % Test changing ()
