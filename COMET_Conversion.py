import pandas as pd

# Define the file path for gene x cell
csv_file_matrix = 'C:/Users/Stevens Lab/Desktop/LanserLIGER/seurat3xScaled_GeneXCellData.csv'

# Load the count matrix as a pandas dataframe
df = pd.read_csv(csv_file_matrix,)
print(df.head())
# # Save the count matrix as a tab-separated .txt file
#df.to_csv( 'C:/Users/Stevens Lab/Desktop/LanserLIGER/COMET/seurat3xScaled_GeneXCellData.txt', header=True, index = True, sep='\t')
#
# #define the path for UMAP coordinates
# csv_file_UMAP = 'C:/Users/Stevens Lab/Desktop/LanserLIGER/UMAP_INMFCoordindates.csv'
#
# # Load the  UMAP coordinates as a pandas dataframe
# count_UMAP = pd.read_csv(csv_file_UMAP, header=0, index_col=0)
#
# # Save the  UMAP coordinates as a tab-separated .txt file
# count_UMAP.to_csv( 'C:/Users/Stevens Lab/Desktop/LanserLIGER/COMET/UMAP_INMFCoordindates.txt', sep='\t')
#
# #define the path for cluster file
# csv_file_cluster = 'C:/Users/Stevens Lab/Desktop/LanserLIGER/seurat3x_active.ident.csv'
#
# # Load the cluster file as a pandas dataframe
# count_cluster = pd.read_csv(csv_file_cluster, header=0, index_col=0)
#
# # Save the cluster fileas a tab-separated .txt file
# count_cluster.to_csv( 'C:/Users/Stevens Lab/Desktop/LanserLIGER/COMET/seurat3x_active.ident.txt', sep='\t')
