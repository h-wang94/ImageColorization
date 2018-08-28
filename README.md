IDEAS:
-Gabor filter parametrizado pelo tamanho da imagem
-nClusters (cor) ajustado pelo erro de recolorização.
-Edges no relabeling
-Clustering: Incluir um canal de textura (emular a classificação manual)

TODO:
-Remover todas as referências à validSuperpixels (nova abordagem é robusta).
-Refactoring:
	-Parei em: Source sampling
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