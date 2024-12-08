#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include<QLabel>
#include <QMainWindow>
#include <QVector>
#include <QLineEdit>
#include <QPushButton>
#include <QDebug>
#include <QString>
#include <QChart>
#include <QChartView>
#include <QLineSeries>
#include <QVBoxLayout>
#include <QValueAxis>
#include <cmath>
#include <QImage>

QT_BEGIN_NAMESPACE


namespace Ui
{
class MainWindow;
}
QT_END_NAMESPACE

// Объявление структуры cmplx для хранения комплексных чисел
typedef struct cmplx {
    float real;
    float image;
} cmplx;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
public slots:

signals:

private slots:

    void onPlayButtonClicked();
    void onParametrsButtonClicked();

private:

    void fourea (cmplx *data, int n, int is);
    void setupChart();


    Ui::MainWindow *ui;
    QChart *ishod;
    QChartView *chartView;
    QLineSeries *series;
    QValueAxis *AxisT1;
    QValueAxis *AxisX1;

    QChart *fft;
    QChartView *fftchartView;
    QLineSeries *fftseries;
    QLineSeries *fftcleanedseries;
    QValueAxis *AxisT2;
    QValueAxis *AxisX2;

    QLineSeries *noiseseries;

    QChart *CleanClear;
    QChartView *CleanClearChartView;
    QLineSeries *CleanClearSeries;
    QValueAxis *AxisT3;
    QValueAxis *AxisX3;

    QLineSeries *cleanfftseries;

    QLabel *l;
    QLabel *l1;
    QLabel *l2;
    QLabel *l3;
    QLabel *l4;
    QLabel *l5;
    QLabel *l6;
    QLabel *l7;
    QLabel *l8;
    QLabel *l9;
    QLabel *l10;
    QLabel *l11;
    QLabel *l12;
    QLabel *l13;

    QLineEdit *bA1;
    QLineEdit *cf1;
    QLineEdit *da1;
    QLineEdit *bA2;
    QLineEdit *cf2;
    QLineEdit *da2;
    QLineEdit *bA3;
    QLineEdit *cf3;
    QLineEdit *da3;
    QLineEdit *fdL;
    QLineEdit *Num;
    QLineEdit *alphaL;
    QLineEdit *gammaL;
    QLineEdit *deltaL;
};
#endif // MAINWINDOW_H
