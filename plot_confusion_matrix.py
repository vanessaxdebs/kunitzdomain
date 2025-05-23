import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Example confusion matrix counts (replace these with your actual results)
TP = 58
FP = 0
FN = 13
TN = 29

# Build the confusion matrix
conf_matrix = np.array([[TN, FP],
                        [FN, TP]])

# Calculate percentages
total = conf_matrix.sum()
percent_matrix = 100 * conf_matrix / total

# Prepare annotation labels (e.g., "29.03%")
labels = np.array([["{:.2f}%".format(percent_matrix[i, j]) for j in range(2)] for i in range(2)])

# Create the 'images' directory if it doesn't exist
if not os.path.exists('images'):
    os.makedirs('images')

# Plot
plt.figure(figsize=(6,5))
sns.heatmap(percent_matrix, annot=labels, fmt='', cmap='Blues', cbar=True, 
            xticklabels=['Pred 0', 'Pred 1'], yticklabels=['True 0', 'True 1'])
plt.xlabel('Predicted label')
plt.ylabel('True label')
plt.title('Confusion Matrix (percentages)')
plt.tight_layout()

# Save the figure to the 'images' directory
plt.savefig('images/confusion_matrix_kunitz.png', dpi=150)

# Since this is a headless environment, no need to call plt.show()
# plt.show()  # This line is commented out or removed

print("Confusion matrix plot saved to 'images/confusion_matrix_kunitz.png'")

