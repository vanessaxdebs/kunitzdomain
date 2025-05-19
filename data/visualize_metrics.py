import matplotlib.pyplot as plt

def plot_confusion_bar(metrics, output_path):
    """Plot confusion matrix as a bar chart."""
    labels = ['TP', 'FP', 'FN', 'TN']
    values = [metrics['TP'], metrics['FP'], metrics['FN'], metrics['TN']]
    plt.figure(figsize=(6,4))
    bars = plt.bar(labels, values, color=['green', 'red', 'orange', 'blue'])
    plt.title("Confusion Matrix")
    plt.ylabel("Count")
    plt.xlabel("Class")
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x()+0.25, yval+0.1, int(yval))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_metrics_summary(metrics, output_path):
    """Plot precision, recall, accuracy, F1 as a bar chart."""
    labels = ['Accuracy', 'Precision', 'Recall', 'F1']
    values = [metrics['accuracy'], metrics['precision'], metrics['recall'], metrics['f1']]
    plt.figure(figsize=(6,4))
    bars = plt.bar(labels, values, color=['skyblue', 'lightgreen', 'gold', 'violet'])
    plt.ylim(0,1)
    plt.title("Performance Metrics")
    plt.ylabel("Score")
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x()+0.15, yval+0.02, f"{yval:.2f}")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
