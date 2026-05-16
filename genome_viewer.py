# Genome View Backend Code
import sys
import bamnostic
from PyQt5.QtWidgets import (
    QApplication, QMainWindow,
    QGraphicsView, QGraphicsScene,
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QSpinBox, QPushButton,
    QGraphicsRectItem, QGraphicsSimpleTextItem
)
from PyQt5.QtGui import QColor, QBrush, QPen, QFont
from PyQt5.QtCore import QRectF, Qt
# GENOME VIEWER
class GenomeViewer(QGraphicsView):
    def __init__(self, fasta_path, bam_path):
        super().__init__()
        self.reference = self.load_fasta(fasta_path)
        self.bam_path = bam_path
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.base_width = 15
        self.row_height = 22
        self.ref_height = 25
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.font = QFont("Courier", 8)
    # FASTA LOADER
    def load_fasta(self, fasta_path):
        seq = []
        with open(fasta_path) as f:
            for line in f:
                if not line.startswith(">"):
                    seq.append(line.strip().upper())
        return "".join(seq)
    # FETCH ALL READS
    def fetch_all_reads(self):
        reads = []
        with bamnostic.AlignmentFile(self.bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                reads.append(read)
        return reads
    # IGV STACKING
    def stack_reads(self, reads, start, end):
        rows_end = []
        layout = []
        for read in reads:
            if read.is_unmapped:
                continue
            rstart = read.reference_start + 1
            rend = read.reference_end  # ✅ FIXED (CIGAR-aware)
            if rend < start or rstart > end:
                continue
            placed = False
            for i in range(len(rows_end)):
                if rstart >= rows_end[i]:
                    rows_end[i] = rend
                    layout.append((i, read))
                    placed = True
                    break
            if not placed:
                rows_end.append(rend)
                layout.append((len(rows_end) - 1, read))
        return layout
    # BASE COLOR
    def get_base_color(self, base):
        colors = {
            "A": QColor(46, 204, 113),
            "T": QColor(231, 76, 60),
            "U": QColor(231, 76, 60),
            "C": QColor(52, 152, 219),
            "G": QColor(243, 156, 18),
            "N": QColor(189, 195, 199)
        }
        return colors.get(base, QColor("gray"))
    # RENDER REGION
    def render_region(self, start, end):
        self.scene.clear()
        reads = self.fetch_all_reads()
        layout = self.stack_reads(reads, start, end)
        # DRAW REFERENCE BASES
        for pos in range(start, min(end, len(self.reference)) + 1):
            base = self.reference[pos - 1]
            x = (pos - start) * self.base_width
            rect = QGraphicsRectItem(
                QRectF(x, 0, self.base_width, self.ref_height)
            )
            rect.setBrush(QBrush(self.get_base_color(base)))
            rect.setPen(QPen(Qt.black))
            self.scene.addItem(rect)
            text = QGraphicsSimpleTextItem(base)
            text.setFont(self.font)
            text.setBrush(QBrush(Qt.white))
            text.setPos(x + 3, 4)
            self.scene.addItem(text)
        # DRAW READS WITH CIGAR
        for row_index, read in layout:
            y = self.ref_height + 10 + row_index * self.row_height
            ref_pointer = read.reference_start + 1
            seq_pointer = 0
            for length, op in read.cigartuples:
                if op == 0:  # M
                    for i in range(length):
                        genomic_pos = ref_pointer + i
                        if genomic_pos < start or genomic_pos >= end:
                            continue
                        base = read.query_sequence[seq_pointer + i]
                        ref_base = self.reference[genomic_pos - 1]
                        x = (genomic_pos - start) * self.base_width
                        rect = QGraphicsRectItem(
                            QRectF(x, y, self.base_width, self.row_height - 4)
                        )
                        # mismatch highlight
                        if base != ref_base:
                            rect.setBrush(QBrush(QColor(155, 89, 182)))  # purple
                        else:
                            rect.setBrush(QBrush(self.get_base_color(base)))
                        rect.setPen(QPen(Qt.black))
                        self.scene.addItem(rect)
                        text = QGraphicsSimpleTextItem(base)
                        text.setFont(self.font)
                        text.setBrush(QBrush(Qt.white))
                        text.setPos(x + 3, y + 3)
                        self.scene.addItem(text)
                    ref_pointer += length
                    seq_pointer += length
                elif op == 1:  # I
                    for i in range(length):
                        genomic_pos = ref_pointer
                        x = (genomic_pos - start) * self.base_width
                        rect = QGraphicsRectItem(
                            QRectF(x, y, self.base_width, self.row_height - 4)
                        )
                        rect.setBrush(QBrush(QColor(142, 68, 173)))  # insertion
                        rect.setPen(QPen(Qt.black))
                        self.scene.addItem(rect)
                    seq_pointer += length
                elif op == 2:  # D
                    for i in range(length):
                        genomic_pos = ref_pointer + i
                        x = (genomic_pos - start) * self.base_width
                        rect = QGraphicsRectItem(
                            QRectF(x, y, self.base_width, self.row_height - 4)
                        )
                        rect.setBrush(QBrush(QColor(0, 0, 0)))  # deletion
                        rect.setPen(QPen(Qt.black))
                        self.scene.addItem(rect)
                    ref_pointer += length
                elif op == 4:  # Soft clip
                    seq_pointer += length
            # Strand arrow
            arrow = QGraphicsSimpleTextItem("←" if read.is_reverse else "→")
            arrow.setFont(self.font)
            arrow.setPos(
                (read.reference_start + 1 - start) * self.base_width,
                y
            )
            self.scene.addItem(arrow)
        # DRAW UNMAPPED READS
        unmapped_row = len(layout) + 2
        for read in reads:
            if read.is_unmapped:
                y = self.ref_height + 10 + unmapped_row * self.row_height
                x = 0
                width = read.query_length * self.base_width
                rect = QGraphicsRectItem(
                    QRectF(x, y, width, self.row_height - 4)
                )
                rect.setBrush(QBrush(QColor(255, 182, 193)))
                rect.setPen(QPen(Qt.black))
                self.scene.addItem(rect)
                unmapped_row += 1
        self.setSceneRect(self.scene.itemsBoundingRect())
    # ZOOM
    def wheelEvent(self, event):
        factor = 1.2
        if event.angleDelta().y() > 0:
            self.scale(factor, 1)
        else:
            self.scale(1 / factor, 1)

# MAIN WINDOW
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VariantXplorer - Advanced IGV Engine")
        self.resize(1200, 700)
        container = QWidget()
        layout = QVBoxLayout()
        region_layout = QHBoxLayout()
        region_layout.addWidget(QLabel("Start:"))
        self.start_box = QSpinBox()
        self.start_box.setMaximum(1000000)
        self.start_box.setValue(100)
        region_layout.addWidget(self.start_box)
        region_layout.addWidget(QLabel("End:"))
        self.end_box = QSpinBox()
        self.end_box.setMaximum(1000000)
        self.end_box.setValue(200)
        region_layout.addWidget(self.end_box)
        self.load_btn = QPushButton("Load Region")
        region_layout.addWidget(self.load_btn)
        layout.addLayout(region_layout)
        self.viewer = GenomeViewer(
            fasta_path="reference.fasta",
            bam_path="aligned_sorted.bam"
        )     
        layout.addWidget(self.viewer)
        container.setLayout(layout)
        self.setCentralWidget(container)
        self.load_btn.clicked.connect(self.load_region)
        self.viewer.render_region(100, 200)

    def load_region(self):
        self.viewer.render_region(
            self.start_box.value(),
            self.end_box.value()
        )
