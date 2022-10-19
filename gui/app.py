from PyQt5.QtWidgets import QWidget, QApplication, QPushButton

application = QApplication([])

mainWindow = QWidget()
mainWindow.setGeometry(0, 0, 350, 400)
mainWindow.setWindowTitle('Button Widget')

pushButton = QPushButton(parent=mainWindow, text='Click me')

mainWindow.show()
application.exec()
