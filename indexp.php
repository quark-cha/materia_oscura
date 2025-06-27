<!DOCTYPE html>
<html>
<head>
  <meta http-equiv="content-type" content="text/html; charset=UTF-8">
  <title>Mis teorias</title>
</head>
<body>
  <p align="center">Mis teorías</p>
  <p align="center">&#x1f10f;©® 2023 Victor Estrada Diaz &#9730;</p><br>
  <div style="width:300px; margin:0 auto;">
    <strong>
      <p></p>
      <div id="col1">
        <p></p>
      </div>
      <div id="col2">
      </div>
    </strong>
  </div>
  <p align="center">
    <?php
    $pdf_dir = "./"; // Directorio actual donde se encuentran los archivos PDF
    $files = glob($pdf_dir . "*.pdf"); // Obtener todos los archivos PDF en el directorio

    foreach ($files as $file) {
      $pdf = basename($file); // Obtener el nombre del archivo con extensión
      echo '<a href="' . $pdf . '" target="_blank">' . $pdf . '</a><br>'; // Mostrar el enlace al archivo PDF
    }
    ?>
  </p>
</body>
</html>
