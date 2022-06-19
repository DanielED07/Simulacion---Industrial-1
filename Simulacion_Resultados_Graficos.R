load("datos_plot.RData")

library(plotly)
library(htmlwidgets)
library(webshot2)

vline <- function(x = 0, color = "black") {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash="dot")
  )
}

plot_ANI <- plot_ly(data=datos_plot,x=~L,y=~ANI,  color = ~Shift) %>% 
  add_lines() %>% 
  layout(shapes = list(vline(x=3)),
         title = "Numero promedio de iteraciones vs L <br> n = 5, contaminación = 4%")

plot_FAP <-plot_ly(data=datos_plot,x=~L,y=~FAP,  color = ~Shift) %>% 
  add_lines() %>% 
  layout(shapes = list(vline(x=3)),
         title = "Porcentaje de alarma falsa vs L <br> n = 5, contaminación = 4%")

plot_TAP <-plot_ly(data=datos_plot,x=~L,y=~TAP,  color = ~Shift) %>% 
  add_lines() %>% 
  layout(shapes = list(vline(x=3)),
         title = "Porcentaje de alarma verdadera vs L <br> n = 5, contaminación = 4%")

plot_MSE_X <-plot_ly(data=datos_plot,x=~L,y=~MSE_X,  color = ~Shift) %>% 
  add_lines() %>% 
  layout(shapes = list(vline(x=3)),
         title ="MSE de las estimaciones de la media vs L <br> n = 5, contaminación = 4%" ,
         yaxis = list(title="MSE"))
  

plot_MSE_Sc4 <-plot_ly(data=datos_plot,x=~L,y=~MSE_Sc4,  color = ~Shift) %>% 
  add_lines() %>% 
  layout(shapes = list(vline(x=3)),
         title ="MSE de las estimaciones de la desviación estandar vs L <br> n = 5, contaminación = 4%" ,
         yaxis = list(title="MSE")) 


saveWidget(widget = plot_ANI, file = "plot_ANI.html")
saveWidget(widget = plot_FAP, file = "plot_FAP.html")
saveWidget(widget = plot_TAP, file = "plot_TAP.html")
saveWidget(widget = plot_MSE_X, file = "plot_MSE_X.html")
saveWidget(widget = plot_MSE_Sc4, file = "plot_MSE_Sc4.html")

webshot(url = "plot_ANI.html", file = "plot_ANI.png", delay = 1, zoom = 4, vheight = 500)
webshot(url = "plot_FAP.html", file = "plot_FAP.png", delay = 1, zoom = 4, vheight = 500)
webshot(url = "plot_TAP.html", file = "plot_TAP.png", delay = 1, zoom = 4, vheight = 500)
webshot(url = "plot_MSE_X.html", file = "plot_MSE_X.png", delay = 1, zoom = 4, vheight = 500)
webshot(url = "plot_MSE_Sc4.html", file = "plot_MSE_Sc4.png", delay = 1, zoom = 4, vheight = 500)

