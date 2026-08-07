

Prototipo de Visor de Imágenes Dicom-CT
===============================
****
|Autor|Zhou Dejun|
|---|---
|Correo electrónico|zhoudejunchn@gmail.com
****
Este es un prototipo de un visor de imágenes Dicom-CT. Las funciones utilizadas en la GUI se describen en dicomrotate.py. Para ver los detalles de las funciones, lea el archivo dicomrotate.py y los comentarios en él. Esta introducción se centra en la GUI, especialmente en la función de los botones. Por favor, contácteme si tiene alguna sugerencia o pregunta.

<img src='https://github.com/McZhD/dicomViewer/blob/master/ScreenShot.png'>

# 1.How to start.
Este programa está basado en Python3, los módulos requeridos se listan en el archivo requirements.txt. Para ejecutar este programa, coloque viewergui.py y dicomanal.py en la misma carpeta, acceda a esta carpeta y ejecute el comando "pythonw viewergui.py". Tenga en cuenta que se utiliza "pythonw" en lugar de "python" para evitar el problema "This program needs access to the screen".

# 2.Open Patient Button
Haga clic en este botón y seleccione la carpeta que contiene los archivos '.dcm' del paciente. A continuación, la ruta de la carpeta se mostrará en el cuadro de texto adyacente y se visualizará la imagen del primer archivo .dcm (ordenado por nombre de archivo).

# 3.Previous Slice Button
Haga clic en este botón para mostrar la imagen anterior del archivo .dcm en la lista de archivos. Si se está mostrando el primer corte, haga clic en este botón y se mostrará el último corte.

# 4.Next Slice Button
Haga clic en este botón para mostrar la imagen siguiente del archivo .dcm en la lista de archivos. Si se está mostrando el último corte, haga clic en este botón y se mostrará el primer corte.

# 5.Go Button
Ingrese el número de corte en el cuadro de texto anterior, luego haga clic en este botón y se mostrará el corte solicitado.

# 6.Show Structure Check Box
Si esta casilla está marcada, la región de la estructura se mostrará en la imagen. La leyenda de la estructura se visualizarnsestrará en el lado izquierdo de la imagen.

# 7.Show Isodose Check Box
Si esta casilla está marcada, la distribución de la dosis se mostrará en la imagen. La leyenda de la isodosis se /nsestrará en el lado derecho de la imagen. Tenga en cuenta que si no hay distribución de dosis en la imagen, no cambiará nada al marcar esta casilla.
