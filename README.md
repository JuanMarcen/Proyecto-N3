# Generador de lluvia horaria

Software desarrollado para la construcción de un generador de precipitación horaria en el área del incendio de Ateca. Esta tarea forma parte del proyecto de investigación HIDROGIF: Desarrollo de una plataforma computacional para el análisis integral del riesgo de
degradación hidrológico-forestal por grandes incendios bajo escenarios de cambio climático.

<img width="1659" height="955" alt="HIDROGIF" src="https://github.com/user-attachments/assets/298002f1-a086-41fe-a229-73cee086af2b" />

## Estructura scripts

El lenguaje de programación utilizado es R. Los scripts se estructuran de la siguiente forma

### Organización de los datos

- Creación de datos horarios y diarios a partir de los datos 15 minutales proporcionados por la Confederación Hidrográfica del Ebro (CHEBRO): [01_dataframes.R](01_dataframes.R)
- Obtención de datos climáticos diarios del reanálisis ERA5: [04_Dataframes_ERA5.R](04_Dataframes_ERA5.R)
- Creación del conjunto de datos de las estaciones a utilizar (ubicación, nombre...): [03_stations.R](03_stations.R)
- Creación de mapas: [mapa_points.R](mapa_points.R)

### Exploratorio

- Análisis exploratorio de los datos: [Exploratorio.R](Exploratorio/Exploratorio.R)

### Modelización

- Exploración de los mejores modelos locales:
  - Modelo de ocurrencia horario (MHO): [02.1_MHO.R](02.1_MHO.R)
  - Modelo de cantidad horario (MHQ): [02.2_MHQ.R](02.2_MHQ.R)
  - Modelo de ocurrencia diario (MDO): [05.1_MDO.R](05.1_MDO.R)
  - Modelo de cantidad diario (MDO): [05.2_MDQ.R](05.2_MDQ.R)
 
- Estudio de los mejores modelos comunes a todas estaciones: [07_modeloComun.R](07_modeloComun.R)
- Análisis de bondad de los modelos según métricas: [06_comp_by_methods.R](06_comp_by_methods.R), [methods.R](methods.R)

### Generadores de lluvia
- Creación de generadores de lluvia y análisis de sensibilidad: [09_An.sensibilidad.R](09_An.sensibilidad.R)
