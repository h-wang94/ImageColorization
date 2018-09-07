TODO:
-Gabor filter parametrizado pelo tamanho da imagem
-nClusters (cor) ajustado pelo erro de recolorização.
-Refactoring:
	-Parei em: Source sampling
	-Nova versão que só tenha modos de operação para superpixels.
-Color Refinement: Ajuste sat do código de Gupta
-Remover estruturas de sampling

Code:
-Características full e salvar.
	-Na extração quando carregado, faz só o subset.
	-save coloca a configuração utilizada.
-Padronizar as estruturas do código 
	-source, target, samples, clusters 
-Explicitar os parâmetros escondidos
	-Quando passo a estrutura inteira para uma função e uso internamente os atributos.
-Indexar todas as chamadas de "figure"