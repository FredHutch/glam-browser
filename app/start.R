library("shiny")
options(shiny.host="0.0.0.0", shiny.port=7777)
runApp(Sys.getenv("APP_DIR"))
