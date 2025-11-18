# Generador de lluvia horaria

Software desarrollado para la construcción de un generador de precipitación horaria en el área del incendio de Ateca. Esta tarea forma parte del proyecto de investigación HIDROGIF: Desarrollo de una plataforma computacional para el análisis integral del riesgo de
degradación hidrológico-forestal por grandes incendios bajo escenarios de cambio climático.

<img width="1659" height="955" alt="HIDROGIF" src="https://github.com/user-attachments/assets/298002f1-a086-41fe-a229-73cee086af2b" />

## Estructura scripts

El lenguaje de programación utilizado es R. Los scripts se estructuran de la siguiente forma

### Organización de los datos

- Creación de datos horarios y diarios a partir de los datos 15 minutales proporcionados por la Confederación Hidrográfica del Ebro (CHEBRO): 01_dataframes.R
- Obtención de datos climáticos diarios del reanálisis ERA5
- Creación del conjunto de datos de las estaciones a utilizar (ubicación, nombre...)

### Exploratorio

- Análisis exploratorio de los datos

### Modelización

- Exploración de los mejores modelos locales:
  - Modelo de ocurrencia horario (MHO):
  - Modelo de cantidad horario (MHQ):
  - Modelo de ocurrencia diario (MDO):
  - Modelo de cantidad diario (MDO):
 
- Estudio de los mejores modelos comunes a todas estaciones:
- Análisis de bondad de los modelos según métricas:

### Generadores de lluvia
- Creación de generadores de lluvia y análisis de sensibilidad:
