#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)

    , ishod(new QChart())
    , chartView(new QChartView(this))
    , series(new QLineSeries())
    , AxisT1(new QValueAxis())
    , AxisX1(new QValueAxis())

    , fft(new QChart())
    , fftchartView(new QChartView(this))
    , fftseries(new QLineSeries())
    , fftcleanedseries(new QLineSeries(this))
    , AxisT2(new QValueAxis())
    , AxisX2(new QValueAxis())

    , noiseseries(new QLineSeries(this))

    , CleanClear(new QChart())
    , CleanClearChartView (new QChartView(this))
    , CleanClearSeries(new QLineSeries())
    , AxisT3(new QValueAxis())
    , AxisX3(new QValueAxis())

    , cleanfftseries(new QLineSeries())
{
    ui->setupUi(this);

    l = new QLabel(this);
    l->setText("A1");
    l->setGeometry(40,30,70,40);

    l1 = new QLabel(this);
    l1->setText("A2");
    l1->setGeometry(160,30,70,40);

    l2 = new QLabel(this);
    l2->setText("A3");
    l2->setGeometry(260,30,70,40);

    l3 = new QLabel(this);
    l3->setText("f1");
    l3->setGeometry(40,90,70,40);

    l4 = new QLabel(this);
    l4->setText("f2");
    l4->setGeometry(160,90,70,40);

    l5 = new QLabel(this);
    l5->setText("f3");
    l5->setGeometry(260,90,70,40);

    l6 = new QLabel(this);
    l6->setText("a1");
    l6->setGeometry(40,150,70,40);

    l7 = new QLabel(this);
    l7->setText("a2");
    l7->setGeometry(160,150,70,40);

    l8 = new QLabel(this);
    l8->setText("a3");
    l8->setGeometry(260,150,70,40);

    l9 = new QLabel(this);
    l9->setText("fd");
    l9->setGeometry(40,210,70,40);

    l10 = new QLabel(this);
    l10->setText("Отсчёты");
    l10->setGeometry(15,270,70,40);

    l11 = new QLabel(this);
    l11-> setText("alpha\n[0,1]");
    l11->setGeometry(175,270,90,50);

    l12 = new QLabel(this);
    l12->setText("gamma\n  [0,1]");
    l12->setGeometry(175,210,90,50);

    l13 = new QLabel(this);
    l13->setText("delta");
    l13->setGeometry(320,240,70,50);

    bA1 = new QLineEdit(this);
    bA1->setGeometry(80,30,70,50);

    cf1 = new QLineEdit(this);
    cf1->setGeometry(80,90,70,50);

    da1 = new QLineEdit(this);
    da1->setGeometry(80,150,70,50);

    bA2 = new QLineEdit(this);
    bA2->setGeometry(180,30,70,50);

    cf2 = new QLineEdit(this);
    cf2->setGeometry(180,90,70,50);

    da2 = new QLineEdit(this);
    da2->setGeometry(180,150,70,50);

    bA3 = new QLineEdit(this);
    bA3->setGeometry(280,30,70,50);

    cf3 = new QLineEdit(this);
    cf3->setGeometry(280,90,70,50);

    da3 = new QLineEdit(this);
    da3->setGeometry(280,150,70,50);

    fdL = new QLineEdit(this);
    fdL->setGeometry(80,210,70,50);

    Num = new QLineEdit(this);
    Num->setGeometry(80,270,70,50);

    alphaL = new QLineEdit(this);
    alphaL-> setGeometry(230,270,70,50);

    gammaL = new QLineEdit(this);
    gammaL->setGeometry(230,210,70,50);

    deltaL = new QLineEdit(this);
    deltaL->setGeometry(360,240,70,50);

    chartView->setGeometry(500,1,800,310);
    chartView->setPalette(Qt::white);
    chartView->setChart(ishod);
    ishod->addSeries(series);

    ishod->setTitle("График исходного сигнала");
    AxisT1 = new QValueAxis();
    AxisX1 = new QValueAxis();
    AxisT1->setTitleText("T, c");
    AxisX1->setTitleText("X");
    ishod->addAxis(AxisT1, Qt::AlignBottom);
    ishod->addAxis(AxisX1, Qt::AlignLeft);
    series->attachAxis(AxisX1);
    series->attachAxis(AxisT1);

    ishod->addSeries(noiseseries);
    noiseseries->attachAxis(AxisX1);
    noiseseries->attachAxis(AxisT1);


    fftchartView->setGeometry(500,281,800,313);
    fftchartView->setPalette(Qt::white);
    fftchartView->setChart(fft);
    fft->addSeries(fftseries);

    fft->setTitle("График аплитудного спектра");
    AxisT2 = new QValueAxis();
    AxisX2 = new QValueAxis();
    AxisT2->setTitleText("f, Герц");
    AxisX2->setTitleText("An");
    fft->addAxis(AxisT2, Qt::AlignBottom);
    fft->addAxis(AxisX2, Qt::AlignLeft);
    fftseries->attachAxis(AxisX2);
    fftseries->attachAxis(AxisT2);


    CleanClearChartView->setGeometry(500,561,800,310);
    CleanClearChartView->setPalette(Qt::white);
    CleanClearChartView->setChart(CleanClear);
    CleanClear->addSeries(CleanClearSeries);

    CleanClear->setTitle("График исходного и очищенного сигнала");
    AxisT3 = new QValueAxis();
    AxisX3 = new QValueAxis();
    AxisT3->setTitleText("Т, с");
    AxisX3->setTitleText("An");
    CleanClear->addAxis(AxisT3, Qt::AlignBottom);
    CleanClear->addAxis(AxisX3, Qt::AlignLeft);
    CleanClearSeries->attachAxis(AxisX3);
    CleanClearSeries->attachAxis(AxisT3);

    QPushButton *ParametrsButton = new QPushButton("По умолчанию",this);
    ParametrsButton->setGeometry(160, 350, 100, 40);
    ParametrsButton->setStyleSheet(
        "QPushButton { background-color: green; color: white; }"
        "QPushButton:hover { background-color: darkgreen; }"
        );
    connect(ParametrsButton, &QPushButton::clicked,this, &MainWindow::onParametrsButtonClicked);

    QPushButton *playButton = new QPushButton("Старт", this);
    playButton->setGeometry(80, 350, 60, 40);
    playButton->setStyleSheet(
        "QPushButton { background-color: blue; color: white; }"
        "QPushButton:hover { background-color: darkblue; }"
        );
    connect(playButton, &QPushButton::clicked, this, &MainWindow::onPlayButtonClicked);

}

MainWindow::~MainWindow()
{
    delete ui;
    delete ishod;
    delete chartView;
    delete series;
    delete fft;
    delete fftchartView;
    delete fftseries;
    delete noiseseries;
    delete CleanClear;
    delete CleanClearSeries;
    delete CleanClearChartView;
    delete fftcleanedseries;
    delete AxisT1;
    delete AxisX1;
    delete AxisX2;
    delete AxisT2;
    delete AxisX3;
    delete AxisT3;

}

void MainWindow::fourea(cmplx *data,int n,int is)
{

    int i,j,istep;
    int m,mmax;
    float r,r1,theta,w_r,w_i,temp_r,temp_i;
    float pi=3.1415926f;

    r=pi*is;
    j=0;
    for(i=0;i<n;i++)
    {
        if(i<j)
        {
            temp_r=data[j].real;
            temp_i=data[j].image;
            data[j].real=data[i].real;
            data[j].image=data[i].image;
            data[i].real=temp_r;
            data[i].image=temp_i;
        }
        m=n>>1;
        while(j>=m)
        {
            j-=m; m=(m+1)/2;
        }
        j+=m;
    }
    mmax=1;
    while(mmax<n)
    {
        istep=mmax<<1;
        r1=r/(float)mmax;
        for(m=0;m<mmax;m++)
        {
            theta=r1*m;
            w_r=(float)cos((double)theta);
            w_i=(float)sin((double)theta);
            for(i=m;i<n;i+=istep)
            {
                j=i+mmax;
                temp_r=w_r*data[j].real - w_i*data[j].image;
                temp_i=w_r*data[j].image + w_i*data[j].real;
                data[j].real=data[i].real - temp_r;
                data[j].image=data[i].image - temp_i;
                data[i].real+=temp_r;
                data[i].image+=temp_i;
            }
        }
        mmax=istep;
    }
    if(is>0)
        for(i=0;i<n;i++)
        {
            data[i].real/=(float)n;
            data[i].image/=(float)n;
        }

}

void MainWindow::onParametrsButtonClicked()
{
    bA1->setText("10");
    bA2->setText("20");
    bA3->setText("30");
    cf1->setText("10");
    cf2->setText("20");
    cf3->setText("30");
    da1->setText("10");
    da2->setText("20");
    da3->setText("30");
    fdL->setText("300");
    Num->setText("128");
    gammaL->setText("0.5");
    alphaL->setText("1");
}

void MainWindow::onPlayButtonClicked()
{

    auto convertToDouble = [](const QString &input) -> double //реализация арифм опер в кнопках
    {
        QString str = input.trimmed(); // Убираем лишние пробелы

        if (str.contains("M_PI"))
        {
            str.replace("M_PI", QString::number(M_PI)); // Заменяем M_PI на его числовое значение
        } else if (str.contains("PI"))
        {
            str.replace("PI", QString::number(M_PI)); // Для тех, кто вводит просто "PI"
        }

        // Обработка выражений с делением и умножением
        if (str.contains("/"))
        {
            QStringList parts = str.split("/"); // Разделяем по символу "/", split разделяет на 2 части по символу()
            if (parts.size() == 2) {
                double chislitel = parts[0].toDouble();
                double znamenatel = parts[1].toDouble();
                return chislitel / znamenatel;
            }
        } else if (str.contains("*"))
        {
            QStringList parts = str.split("*"); // Разделяем по символу "*"
            if (parts.size() == 2) {
                double firstPart = parts[0].toDouble();
                double secondPart = parts[1].toDouble();
                return firstPart * secondPart;
            }
        }

        return str.toDouble();

    };


    QString A1Str = bA1->text();
    double A1 = A1Str.toDouble();

    QString f1Str = cf1->text();
   double f1 = convertToDouble(f1Str);

    QString a1Str = da1->text();
    double a1 = convertToDouble(a1Str);

    QString A2Str = bA2->text();
    double A2 = A2Str.toDouble();

    QString f2Str = cf2->text();
    double f2 = convertToDouble(f2Str);

    QString a2Str = da2->text();
    double a2 = convertToDouble(a2Str);

    QString A3Str = bA3->text();
    double A3 = A3Str.toDouble();

    QString f3Str = cf3->text();
    double f3 = convertToDouble(f3Str);

    QString a3Str = da3->text();
    double a3 = convertToDouble(a3Str);

    QString fdStr = fdL->text();
    double fd = fdStr.toDouble();

    QString NumStr = Num->text();
    double N = NumStr.toDouble();

    QString alphaStr = alphaL->text();
    double alpha = alphaStr.toDouble();

    QString gammaStr = gammaL->text();
    double gamma = gammaStr.toDouble();

    QString deltaStr = deltaL->text();
    double delta = deltaStr.toDouble();


    double minT = std::numeric_limits<double>::max();
    double maxT = std::numeric_limits<double>::lowest(); // Максимальное значение для оси T (Максимальное минимальное с учетом знака. (отличие от мин)
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();


 qDebug() << "fd значение:" << fd;

    srand(time(0));


    QVector<double> clean;
    QVector<double> noise;
    QVector<double> noisebetta;
    double t;
    clean.resize(N);
    noise.resize(N);
    noisebetta.resize(N);
    series->clear();
    noiseseries->clear();

    double signalSum = 0;
    double noiseValue2=0;
    for(int i=0;i<N;i++)
    {
        t=i/fd;
        clean[i]=A1*cos(2*M_PI*f1*t+a1)+A2*cos(2*M_PI*f2*t+a2)+A3*cos(2*M_PI*f3*t+a3);

        signalSum += pow(clean[i],2); // отсчёты сигнала в квадрате



    double noiseSum=0;
    for(int j=0;j<12;j++)
        {
            double Random = (double) (rand() - RAND_MAX/2); //- рандмакс пополам
            noiseSum+=Random;
        }
    noise[i]=noiseSum; //элемент со случайным значением [n]
    noiseValue2 += pow(noise[i],2);
    }

    double betta = sqrt(alpha*signalSum/noiseValue2);

    for(int i=0;i<N;i++)
    {
        t=i/fd;
        noisebetta[i]=clean[i]+betta*noise[i]; // Сигнал с шумом

        series->append(t,clean[i]);
        noiseseries->append(t, noisebetta[i]);


        if (t < minT) minT = t;
        if (t > maxT) maxT = t;

        if (clean[i] < minX) minX = clean[i];
        if (clean[i] > maxX) maxX = clean[i];

        if (noisebetta[i] < minX) minX = noisebetta[i];
        if (noisebetta[i] > maxX) maxX = noisebetta[i];
    }


    for(int i = 0; i < clean.size(); i++)
    {
        qDebug() << "clean[" << i << "] = " << clean[i];
    }


    for(int i = 0; i < noise.size(); i++)
    {
        qDebug() << "noise[" << i << "] = " << noisebetta[i];
    }


    if (ishod->series().isEmpty())
    {
        ishod->addSeries(series);
        series->attachAxis(AxisT1); // Привязываем ось времени (T)
        series->attachAxis(AxisX1); // Привязываем ось амплитуды (X)
    }



    AxisT1->setRange(minT, maxT); // Устанавливаем диапазон оси T
    AxisX1->setRange(minX, maxX); // Устанавливаем диапазон оси X

    ishod->update();

    //FFT/////////////////////////////////////////////////////////////////////////
    qDebug("BPF");
    //gamma = 1/(1+alpha); //Доля энергии которую xотим оставить

    //delta = разность отсчётов в квадрате, это энергия, делить на сум исх сигнала



    QVector<cmplx> fftData(N);
    for (int i = 0; i < N; i++)
    {
        fftData[i].real = noisebetta[i];
        fftData[i].image = 0;  // Исходные данные реальные, мнимая часть нулевая
    }


    // Применяем БПФ
    fourea(fftData.data(), N,  -1);  // Вызов БПФ

    fftseries->clear();  // Очищаем серию перед построением спектра

    double minFf = 0;
    double maxFf = fd;  // Верхняя граница частотного диапазона
    double minAmpl = std::numeric_limits<double>::max();
    double maxAmpl = std::numeric_limits<double>::lowest();

    // Построение амплитудного спектра

    for (int i = 0; i < N ; i++)
    {

        double Ff = i * fd / N;  // Частота
        double amplitude = sqrt(fftData[i].real * fftData[i].real + fftData[i].image * fftData[i].image);  // Амплитуда

        fftseries->append(Ff, amplitude);

        if (amplitude < minAmpl) minAmpl = amplitude;
        if (amplitude > maxAmpl) maxAmpl = amplitude;
    }


    // Если график амплитудного спектра не был добавлен ранее, добавляем его
    if (fft->series().isEmpty())
    {
        fft->addSeries(fftseries);
        fftseries->attachAxis(AxisT2); // Привязываем ось частоты (T)
        fftseries->attachAxis(AxisX2); // Привязываем ось амплитуды (X)
    }


    // Пороговое значение энергии


    QVector<cmplx> cleanedfftData(N);

    for (int i = 0; i < N; i++)
    {
        cleanedfftData[i].real = noisebetta[i];
        cleanedfftData[i].image = 0;  // Исходные данные реальные, мнимая часть нулевая
    }

    fourea(cleanedfftData.data(), N,  -1);



    double Esum = 0;
    for (int i = 0; i < N; i++)
    {
        Esum += cleanedfftData[i].real * cleanedfftData[i].real + cleanedfftData[i].image * cleanedfftData[i].image;
         // Сумма квадратов значений зашумленного сигнала
    }
    double Es = gamma * Esum; //доля энергии оставляемой гамма,а мне надо  сравнивать с клин пов 2 (энергия чистого сигнала)


    double  Et=0;
    for (int i = 0; i < N; i++)
    {
        Et+=pow(cleanedfftData[i].real,2)+pow(cleanedfftData[N-i -1].real,2)+ pow(cleanedfftData[i].image,2) + pow(cleanedfftData[N-i -1].image,2) ;
        if (Et > Es)
        {
            // Зануляем от i до N-i
            for (int j = i+1; j < N - i-1; j++) //i = th
            {
                cleanedfftData[j].real = 0;
                cleanedfftData[j].image = 0;
            }
            break; // Прекращаем цикл после очистки
        }
    }

    fftcleanedseries->clear();


    double minCleanedAmpl = std::numeric_limits<double>::max();
    double maxCleanedAmpl = std::numeric_limits<double>::lowest();
    for (int i = 0; i < N; i++)
    {
        double Ff = i * fd / N;  // Частота
        double cleanedAmplitude = sqrt(cleanedfftData[i].real * cleanedfftData[i].real + cleanedfftData[i].image * cleanedfftData[i].image);  // Амплитуда очищенного спектра

        fftcleanedseries->append(Ff, cleanedAmplitude);

        if (cleanedAmplitude < minCleanedAmpl) minCleanedAmpl = cleanedAmplitude;
        if (cleanedAmplitude > maxCleanedAmpl) maxCleanedAmpl = cleanedAmplitude;
    }



    // Если график очищенного спектра не был добавлен ранее, добавляем его
    if (fft->series().size() < 2)
    {
        fft->addSeries(fftcleanedseries);
        fftcleanedseries->attachAxis(AxisT2);  // Привязываем ось частоты (T)
        fftcleanedseries->attachAxis(AxisX2);  // Привязываем ось амплитуды (X)
    }

    // Настройка осей для обоих спектров
    AxisT2->setRange(minFf, maxFf);  // Диапазон частот
    AxisX2->setRange(std::min(minAmpl, minCleanedAmpl), std::max(maxAmpl, maxCleanedAmpl));  // Диапазон амплитуд


    fft->update();  // Обновление графика
    qDebug("end BPF");


// ////////CleanClearSignal and cleanfft//////////////////////////////
// исходный сигнал
    QVector<double> CleanClearSignal;
    CleanClearSignal.resize(N);
    CleanClearSeries->clear();

    for(int i=0;i<N;i++)
    {
        t=i/fd;

        CleanClearSeries->append(t,clean[i]);

        if (t < minT) minT = t;
        if (t > maxT) maxT = t;

        if (clean[i] < minX) minX = clean[i];
        if (clean[i] > maxX) maxX = clean[i];

    }

    // очищенный сигнал после обратного фурье преобразования
    QVector<double> cleanfft;
    cleanfft.resize(N);
    cleanfftseries->clear();

    fourea(cleanedfftData.data(), N,  1);

    for (int i = 0; i < N; i++)
    {
        double t = i/fd;
    for (int j = i; j < 2*N - i; j++)
        {


        cleanfft[i] = cleanedfftData[i].real;

        cleanfftseries->append(t, cleanfft[i]);
        }

        if (cleanfft[i] < minX) minX = cleanfft[i];
        if (cleanfft[i] > maxX) maxX = cleanfft[i];
    }

        // Исходный сигнал
    if (CleanClear->series().isEmpty())
    {
        CleanClear->addSeries(CleanClearSeries);

        // Привязываем оси
        CleanClearSeries->attachAxis(AxisT3);
        CleanClearSeries->attachAxis(AxisX3);


    }

        // Сигнал после обратного Фурье
    if (CleanClear->series().size()<2)
    {
    CleanClear->addSeries(cleanfftseries);

        cleanfftseries->attachAxis(AxisT3);
        cleanfftseries->attachAxis(AxisX3);
    }

    // Обновляем диапазоны осей
    AxisT3->setRange(minT, maxT);
    AxisX3->setRange(minX, maxX);


    CleanClear->update();

    double raznostsignalov_2 = 0; //энергия разности сигналов (х-хочищ)
    for (int i=0;i<N;i++)
    {
        raznostsignalov_2 += pow((clean[i]) - (cleanfft[i]),2);
    }

    delta = raznostsignalov_2/signalSum; //делить на энергию исходного сигнала (чистый)

    deltaL->setText(QString::number(delta));

// ////////cleanfft////////////////////////////// подсчитать энергию чистого сигнала (сумма всех отсчётов в квадрате) энергия отклонения до фильтрации и после нее, долен стать меньше после удаления

}
