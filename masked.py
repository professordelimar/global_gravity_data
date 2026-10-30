import csv
from pathlib import Path

# ==========================================
# CONFIGURACOES
# ==========================================

ARQUIVO_ENTRADA = "gravity_results.txt"
ARQUIVO_SAIDA = "output.txt"

# Informe os limites da área
LONG_MIN = -52.1
LONG_MAX = -51.9
LAT_MIN = -3.5
LAT_MAX = -3.2

DELIMITADOR = ","  # None = espaços ou tabulações


# ==========================================
# FUNÇÕES
# ==========================================

def separar_linha(linha):
    if DELIMITADOR is None:
        return linha.strip().split()
    return next(csv.reader([linha], delimiter=DELIMITADOR))


def coordenada_mais_proxima(arquivo, coluna, alvo):
    """
    Procura a coordenada existente mais próxima
    do limite informado, sem carregar o arquivo
    inteiro na memória.
    """
    melhor = None
    menor_distancia = float("inf")

    with open(arquivo, "r", encoding="utf-8-sig") as f:
        cabecalho = separar_linha(next(f))
        indice = cabecalho.index(coluna)

        for linha in f:
            if not linha.strip():
                continue

            campos = separar_linha(linha)

            try:
                valor = float(campos[indice])
            except (ValueError, IndexError):
                continue

            distancia = abs(valor - alvo)

            if distancia < menor_distancia:
                menor_distancia = distancia
                melhor = valor

    return melhor


def extrair_area(arquivo, saida,
                 long_min, long_max,
                 lat_min, lat_max):

    if long_min > long_max or lat_min > lat_max:
        raise ValueError("Os limites mínimos devem ser menores que os máximos.")

    print("Localizando coordenadas mais próximas...")

    lon_min_real = coordenada_mais_proxima(
        arquivo, "longitude", long_min
    )
    lon_max_real = coordenada_mais_proxima(
        arquivo, "longitude", long_max
    )
    lat_min_real = coordenada_mais_proxima(
        arquivo, "latitude", lat_min
    )
    lat_max_real = coordenada_mais_proxima(
        arquivo, "latitude", lat_max
    )

    if None in (lon_min_real, lon_max_real,
                lat_min_real, lat_max_real):
        raise ValueError("Não foram encontradas coordenadas válidas.")

    print("\nLimites utilizados:")
    print(f"Longitude: {lon_min_real} até {lon_max_real}")
    print(f"Latitude:  {lat_min_real} até {lat_max_real}")

    total = 0

    with open(arquivo, "r", encoding="utf-8-sig") as entrada, \
         open(saida, "w", encoding="utf-8") as output:

        cabecalho = separar_linha(next(entrada))

        indice_lon = cabecalho.index("longitude")
        indice_lat = cabecalho.index("latitude")
        indice_excluir = cabecalho.index("orto_grav")

        # Preserva todas as colunas, exceto orto_grav
        indices = [
            i for i in range(len(cabecalho))
            if i != indice_excluir
        ]

        separador_saida = (
            DELIMITADOR if DELIMITADOR is not None else "\t"
        )

        output.write(
            separador_saida.join(cabecalho[i] for i in indices) + "\n"
        )

        print("\nExtraindo dados...")

        for linha in entrada:
            if not linha.strip():
                continue

            campos = separar_linha(linha)

            try:
                lon = float(campos[indice_lon])
                lat = float(campos[indice_lat])
            except (ValueError, IndexError):
                continue

            if (lon_min_real <= lon <= lon_max_real and
                lat_min_real <= lat <= lat_max_real):

                output.write(
                    separador_saida.join(
                        campos[i] for i in indices
                    ) + "\n"
                )

                total += 1

                if total % 100000 == 0:
                    print(f"{total:,} linhas extraídas...")

    print("\nProcessamento concluído!")
    print(f"Total de linhas de dados: {total:,}")
    print(f"Arquivo gerado: {Path(saida).resolve()}")

if __name__ == "__main__":

    extrair_area(
        ARQUIVO_ENTRADA,
        ARQUIVO_SAIDA,
        LONG_MIN,
        LONG_MAX,
        LAT_MIN,
        LAT_MAX
    )