# Homogeneização de dados meteorológicos e cálculo de normas climatológicas padrão
Homogeneização de algumas séries de precipitação e temperaturas extremas de Angola e cálculo das suas Normais Climatológicas Padrão 1991-2020

#### Escrito por: Jose A. Guijarro jaguijarrop@gmail.com em 2023-10-31

##  Introdução
Este breve exemplo de homogeneização e geração das Normais Climatológicas Padrão 1991-2020 foi efectuado como contribuição para o Instituto Nacional de Meteorologia e Geofísica de Angola (INAMET) a pedido de António Manuel Lameira Gaspar, meteorologista daquela instituição.

## Dados

As séries de precipitação diária e temperaturas extremas são fornecidas no formato RClimDex para as cinco estações incluídas no ficheiro de texto separado por tabulação.

Os ficheiros de dados têm nomes compostos pelo código da estação e pela extensão txt. Para aplicar o pacote climatol R
para controlo de qualidade, homogeneização e preenchimento de dados destas séries diárias, temos de converter estes ficheiros para o formato de entrada
formato de entrada climatol. A função rclimdex2climatol será a ferramenta correcta para o efeito, mas a primeira coluna
do ficheiro das estações deve conter o nome completo dos ficheiros de dados.

Depois de criámos os ficheiros de entrada para o climatol, podemos proceder à homogeneização, variável a variável.

## Cálculo das Normais Climatológicas Padrão 1991-2020
Para utilizar a ferramenta CLINO, temos de copiar alguns dos seus ficheiros para o nosso diretório de trabalho atual. 
Depois de descomprimir o arquivo CLINO_tool.tgz, os ficheiros a copiar são CLINO.R e os quatro CLINO_*.csv. 
Dois dos campos csv devem ser editados para instruir o software sobre quais são as variáveis a processar e os nomes dos ficheiros de dados homogeneizados.

Agora estamos prontos para carregar a função CLINO na memória do R e executá-la para calcular as Normais Climatológicas Padrão para 1991-2020 e escrevê-las em ficheiros CSV com o formato exigido pela OMM:

      source('CLINO.R')
      CLINO()
