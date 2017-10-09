#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QTextStream>
#include <QLabel>
#include <QString>
#include <QDebug>
#include "visualization.h"
#include "ui_visualization.h"
#include "points.h"
#include "Polynomials.h"
#include "Sylvester.h"
#include "GenericMatrix.h"
#include "GenericMatrix.cpp"
#include "Solution.h"
#include "Calculations.h"
#include "Multiply.h"

# define BUFFSIZE 200


using namespace std;

int counter=0;

Visualization::Visualization(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Visualization)
{
    ui->setupUi(this);
    widget = new QtGnuplotWidget();
    widget->installEventFilter(this);
    widget->setStatusLabelActive(true);
    instance = new QtGnuplotInstance();
    instance->setWidget(widget);
    solve = new QPushButton("SolveB",this);
    solve->setEnabled(false);
    solve->setGeometry(QRect(QPoint(50, 200),QSize(100, 25)));
    points=-1;
    d=-1;
    point=NULL;
    this->solved=false;
    connect(solve,SIGNAL(clicked()),this,SLOT(solve_clicked()));
}

Visualization::~Visualization()
{
    delete ui;
    delete instance;
    delete widget;
}

void Visualization::print_results(GenericMatrix<double *> * syl, SylvPol<GenericMatrix<double> > *sylp, double  k, char  hidden){

    //print results (output , roots)

    ui->outputTxt->setText("");
    int height = syl->GetHigh();
    QString nui =  "Results  : ";
    nui.append("\n\nThe hidden variable is "+hidden);
    nui.append("\nthe rank of y is "+QString::number(height-1));
    nui.append("\n\nPrinting Sylvester Matrix \n\n");
    for(int i=0;i<syl->GetRows();i++){
        for(int j=0;j<syl->GetColumns();j++){
            for (int k=0;k<syl->GetHigh();k++){
                if(syl->GetMatrix(i,j,k)==0)
                    nui.append("0");
                else{
                    nui.append(QString::number(syl->GetMatrix(i,j,k)));
                    nui.append("y^"+QString::number(k));
                }
                if(k!=height-1)
                    nui.append("+");
             }
             nui.append(" ");
        }
        nui.append("\n");
    }
    nui.append("\n");
    nui.append("\nPrinting Md Matrix \n");
    for (int i=0;i<sylp->GetMatrixI(0)->GetRows();i++){
        for(int j=0;j<sylp->GetMatrixI(0)->GetColumns();j++){
            nui.append(QString::number(sylp->GetMatrixI(height-1)->GetMatrix(i,j)));
            nui.append(" ");
        }
        nui.append("\n");
    }
    ui->outputTxt->setText("");
    if(k<pow(10,7)){								//print and return according to K
        nui.append("\nk = " +QString::number(k));
        nui.append(" Bound: non-singular  Md, standard eigenproblem");
    }else{
        nui.append("\nk = " +QString::number(k));
        nui.append(" Bound: ill-conditioned Md, generalized eigenproblem");
    }
    //qDebug()<<nui;
    ui->outputTxt->setText(nui);
    QFile file("roots.txt");
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
        return;
    }
    QTextStream in(&file);
    QString str4 = in.readAll();
    //qDebug()<<str4;
    ui->rootsTxt->setText(str4);


}

void Visualization::on_actionExit_triggered()
{
    qApp->quit();
}

void Visualization::on_actionOpen_triggered() //open file
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), QString(),
                tr("All Files (*.*)"));

    if (!fileName.isEmpty()) {
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly)) {
            QMessageBox::critical(this, tr("Error"), tr("Could not open file"));
            return;
        }
        QTextStream in(&file);
        QString nui=ui->equationsTxt->toPlainText();
        if(!(nui.isEmpty()))
            nui.append("\n");
        QString str = in.readAll();
        nui.append(str);
        ui->equationsTxt->setText(nui);
        file.close();
    }
}

bool Visualization::eventFilter(QObject *obj, QEvent *event) // clicks(user must click appropriate times message is displayed)
{
    QString Optxt;
    if (event->type() == QEvent::MouseButtonPress)
    {
        if (obj == this->widget) {
            QMouseEvent *mouseEvent = static_cast<QMouseEvent *>(event);
            if ((mouseEvent->button() == Qt::LeftButton )&& (points!=-1)) {
                //*instance<< replot\<echo ;
                Optxt=this->widget->getStatusLabel()->text();
                ui->outputTxt->setText(Optxt);
                QStringList ps = Optxt.split(",",QString::SkipEmptyParts);
                this->point[counter][0]=ps[0].toDouble();
                this->point[counter][1]=ps[1].toDouble();
               // qDebug()<<point[counter][0]<<point[counter][1]<<points<<counter;
                counter++;
               emit replot_points();
                if(counter==this->points){
                    QMessageBox::information(this,tr("Points"),tr("You have seleced k points"));
                    emit pointsClosed();
                    this->widget->close();
                }
            }
        }
    }
    return QObject::eventFilter(obj, event);
}

void Visualization::replot_points(){ //after clicking a point->display it
    QString str1="replot\"<echo '"+QString::number(point[counter-1][0]);
    str1.append(" ");
    QString str2 = str1+QString::number(point[counter-1][1]);
    str2.append("'\"notitle\n");
    //qDebug()<<str2;
    *instance << str2;

}

void Visualization::pointsClosed(){//user has inserted appropriate number of points->creation of polynomial through interpolation
    QString pol1 ;
    VectorXd c;
    int xpow=0;
    int ypow=0;
    int maxpow=0;
    MatrixXd interp;
    interp=createInterpolation(point,d,points);
    //qDebug()<<interp(0,0);
    c = computeKernel(interp);
    //cout << c << endl;
    double * coeffs =new double[c.rows()];
    for(int i=0;i<c.rows();i++)
        coeffs[i]=c(i);
    for(int i =0;i<c.rows();i++){
        if(coeffs[i]==0){
            xpow--;
            ypow++;
        }else{
            if(coeffs[i]>0)
                pol1.append("+");
            pol1.append(QString::number(coeffs[i]));
            if(xpow!=0){
                pol1.append("*x^");
                pol1.append(QString::number(xpow));
            }
            if(ypow!=0){
                pol1.append("*y^");
                pol1.append(QString::number(ypow));
            }
        }
        if(ypow==maxpow){
            maxpow++;
            xpow=maxpow;
            ypow=0;
        }else{
            xpow--;
            ypow++;
        }
    }
    QString nui =ui->equationsTxt->toPlainText(); //display new polynomial
    //qDebug() << nui ;
    if (!nui.isEmpty())
        nui.append("\n");
    nui.append(pol1);
    //qDebug() << nui ;
    ui->equationsTxt->setText(nui);
    //ui->equationsTxt->setText(pol1);
    for(int i=0;i<points;i++)
           delete point[i];
    delete point;
    points=-1;
    d=-1;
    counter=0;
}

void Visualization::on_pushButton_clicked() // if problem solved display plots
{
    if(this->solved){
       //only functions + way to print the points

        /* int x=1;
        widget->show();
        widget->resize(QSize(800,600));
        QString str= ui->equationsTxt->toPlainText();
        QStringList strs = str.split("\n",QString::SkipEmptyParts);
        QString str1=strs[0];
        QString str2=strs[1];
        str1.replace("^","**");
        str2.replace("^","**");
        QString str3="set yrange [-10:10]\nset xrange [-10:10]\nset isosamples 500,500\nf(x,y)=";
        str3.append(str1);
        str3.append("\ng(x,y)=");
        str3.append(str2);
        str3.append("\nset zeroaxis\nset contour\nset cntrparam levels discrete 0\nunset key\nset view 0,0\nunset ztics\nunset surface\nsplot f(x,y) , g(x,y) , ");//einai swsto gia tis funcs
        str3.append("\"roots.txt\" with points nocontour pt 7\n");
        *instance<<str3;
        //qDebug()<<str3;*/


        //multiplot functions-roots
        widget->show();
        widget->resize(QSize(800,600));
        QString str= ui->equationsTxt->toPlainText();
        QStringList strs = str.split("\n",QString::SkipEmptyParts);
        QString str1=strs[0];
        QString str2=strs[1];
        str1.replace("^","**");
        str2.replace("^","**");
        QString str3 = "set size 1,1 \nset multiplot\nset size 0.5,0.5\nset origin 0,0.5\nplot \"roots.txt\" with points pt 7\n";
        str3.append("set origin 0.5,0.5\n");
        str3.append("set yrange [-10:10]\nset xrange [-10:10]\nset isosamples 500,500\nf(x,y)=");
        str3.append(str1);
        str3.append("\ng(x,y)=");
        str3.append(str2);
        str3.append("\nset zeroaxis\nset contour\nset cntrparam levels discrete 0\n\nset view 0,0\nunset ztics\nunset surface\nsplot f(x,y) , g(x,y)\nunset multiplot\n");
        //QFile file("ex.txt");
        //QTextStream in(&file);
        //QString str4 =in.readAll();
        *instance<<str3;
        //qDebug()<<str3;

    }else{
        QMessageBox::information(this,tr("Points"),tr("The problem is not solved yet"));
    }

}

void Visualization::on_actionInsert_triggered() //insert polynomila through points
{
    bool ok;
    QString text=QInputDialog::getText(this, tr("Insert x range"),tr("X range:"), QLineEdit::Normal, "2", &ok);
    if(ok){
        int xrange=text.toInt();
        text=QInputDialog::getText(this, tr("Insert y range"),tr("Y range:"), QLineEdit::Normal, "2", &ok);
        if(ok){
            int yrange=text.toInt();
             text = QInputDialog::getText(this, tr("Insert d"),tr("Max D :"), QLineEdit::Normal, "2", &ok);
             if (ok){
                 int md = text.toInt();
                 this->d=md;
                 this->points=((md+1)*(md+2)/2) - 1;
                 point = new double*[points];
                 for (int i=0;i<points;i++)
                      point[i]= new double[2];

                 widget->show();
                 widget->resize(QSize(800,600));
                 QString str1="set tics scale 0.75\nset xtics 1\nset ytics 1\nset yrange ["+QString::number(-yrange);
                 str1.append(":"+QString::number(yrange));
                 str1.append("]\nset xrange ["+QString::number(-xrange));
                 str1.append(":"+QString::number(xrange));
                 str1.append("]\nset xlabel 'x'\nset ylabel 'y'\nset zeroaxis\nset grid\nplot \"<echo ' '\" notitle\n");
                 *instance << str1;
             }

         }
    }
}

void Visualization::solve_clicked(){ //solve button->solve the problem

   int d1,d2,d3,d4,x,y,height,dmax;
   char hidden;
   double * write;

   this->solved = true;

   QString Eqs = ui->equationsTxt->toPlainText();
   QStringList pols = Eqs.split("\n",QString::SkipEmptyParts);

   GenericMatrix<double> * polm1;
   GenericMatrix<double> * polm2;

   QString pols1=pols[0];
   QString pols2=pols[1];

   char * tmp1= (char *)malloc(sizeof(char)*BUFFSIZE);
   char * tmp2= (char *)malloc(sizeof(char)*BUFFSIZE);

   strcpy(tmp1,pols1.toLocal8Bit().constData());
   tmp1[pols1.toLocal8Bit().length()]='\n';
   findmaxpowers(tmp1,&x,&y);
   d2=x;
   d1=y;
   //qDebug()<<tmp1;
   //qDebug()<<d1<<d2;

   polm1=new GenericMatrix<double> (d1+1,d2+1);
   createpolmatrix(tmp1,tmp1,polm1);

   strcpy(tmp2,pols2.toLocal8Bit().constData());
   tmp2[pols2.toLocal8Bit().length()]='\n';
   findmaxpowers(tmp2,&x,&y);
   d4=x;
   d3=y;
   //qDebug()<<tmp2;
   //qDebug()<<d3<<d4;

   polm2=new GenericMatrix<double> (d3+1,d4+1);
   createpolmatrix(tmp2,tmp2,polm2);

   free(tmp1);
   free(tmp2);

   //polm1->PrintMatrix();
   //polm2->PrintMatrix();
   cout << d1 << d2 << d3 << d4<<endl;
   GenericMatrix<double *> * syl;
   if((d2>=d1&&d4>=d3)||(d2>=d1)&&(d2>=d3)||(d4>=d1)&&(d4>=d3))	//If the rank of x is greater or equal to the rank of y
       {
           height=max(d1,d3)+1;			//Height will be the max rank of y
           dmax=d2+d4;
           syl=new GenericMatrix<double *>(dmax,dmax,height);
           hidden='y';
       }
       else															//The rank of y is greater than the rank of x
       {
           height=max(d4,d2)+1;			//Height will be the max rank of x
           dmax=d1+d3;
           syl=new GenericMatrix<double *>(dmax,dmax,height);
           hidden='x';
       }
       cout << "The hidden variable is " << hidden << endl;
       write=new double[height];				//Array that will be written in Sylvester
       for(int i=0; i<height; i++)				//Sylvester is all zeroes
            write[i]=0;
       for(int i=0; i<dmax; i++)
       {
            for(int j=0; j<dmax; j++)
            {
                syl->SetMatrix(i,j,write);
            }
       }
       createSylvester(polm1,polm2,syl);

       cout << "The hidden variable is " << hidden << endl;
       cout << "There are " << d1 << " rows and " << d2 << " columns in the first matrix" << endl;
       cout << "There are " << d3 << " rows and " << d4 << " columns in the second matrix" << endl;
       height=syl->GetHigh();
       dmax=syl->GetRows();
       cout << "There are " << dmax << " rows and " << dmax << " columns in the Sylvester matrix " ;
       cout << "and the rank of y is " << height-1 << endl;
       cout << "Printing Sylvester" << endl;
       syl->PrintMatrix();
       SylvPol<GenericMatrix <double> > * sylp=new SylvPol<GenericMatrix <double> >(height);
       createPolynomials(syl,sylp);
       int B=7;		//Variables used in the second part
       double bestKfactor,Kfactor;
       int c=0;
       int ts[5];
       ts[0]=rand()%10;
       ts[1]=rand()%10;
       ts[2]=rand()%10;
       ts[3]=rand()%10;
       ts[4]=0;
       //compute K
       if(computeK(sylp,B,height,&Kfactor))				//Calculate K and solve accordingly
       {
           StandardProblem(sylp,dmax,height-1,polm1,polm2,ts,hidden);
       }
       else
       {
           GeneralizedProblem(sylp,dmax,height-1,polm1,polm2,ts,hidden);
       }
       print_results(syl,sylp,Kfactor,hidden);
       bestKfactor=Kfactor;
       //SylvPol<GenericMatrix <double> > * sylpZ=new SylvPol<GenericMatrix <double> >(height);	//Z polynomial
       ts[4]=1;
       for(int tries=0;tries<1;tries++)
       {
           SylvPol<GenericMatrix <double> > * sylpZ=new SylvPol<GenericMatrix <double> >(height);	//Z polynomial
           ts[0]=rand()%10;
           ts[1]=rand()%10;
           ts[2]=rand()%10;
           ts[3]=rand()%10;
           CreateZpol(sylp,sylpZ,ts[0],ts[1],ts[2],ts[3]);	//Create Z polynomial with random ts y=(t1*z+t2)/(t3*z+t4)
           if(computeK(sylpZ,B,height,&Kfactor))		//Compute new K and solve accordingly if K(z)<K(y)
           {
               if(Kfactor<bestKfactor){		//K(z)<K(y)
                   cout<<"yep"<<endl;
                   StandardProblem(sylpZ,dmax,height-1,polm1,polm2,ts,hidden);
                   //break;
               }
                   cout<<"nope"<<endl;		//K(z)>K(y)
          }
          else
          {
               if(Kfactor<bestKfactor){		//K(z)<K(y)
                   GeneralizedProblem(sylpZ,dmax,height-1,polm1,polm2,ts,hidden);
                   cout<<"yep"<<endl;
                   //break;
               }
                   cout<<"nope"<<endl;		//K(z)>K(y)
           }
       delete sylpZ;
       }
       delete polm1;
       delete polm2;
       //delete sylp;
       delete syl;
       delete [] write;
       return ;
}

void Visualization::on_equationsTxt_textChanged()//check if there are two pol/als to solve the problem
{
    QString Eqs = ui->equationsTxt->toPlainText();
    QStringList pols = Eqs.split("\n",QString::SkipEmptyParts);
    //qDebug()<<pols.size();
    if (pols.size()==2){
        this->solve->setEnabled(true);
    }else{
        this->solve->setEnabled(false);
    }
}

void Visualization::on_pushButton_2_clicked()//remove a polynomial
{
    bool ok;
    QString text = QInputDialog::getText(this, tr("Which pol/al do you wish to remove"),tr("No P:"), QLineEdit::Normal, "1", &ok);
    int md = text.toInt();
    QString Eqs = ui->equationsTxt->toPlainText();
    QStringList pols = Eqs.split("\n",QString::SkipEmptyParts);
    QString nui="";
    if((md>pols.size())||(md==0)){
        QMessageBox::information(this,tr("Points"),tr("You selected an invalid pol/mial"));
        return;
    }else{
        for(int i=0;i<pols.size();i++){
            //qDebug()<<nui;
            if(i!=(md-1)){
                if(!(nui.isEmpty()))//kai to teleytaio??\n
                    nui.append("\n");
                nui.append(pols[i]);
            }
        }
    }
    ui->equationsTxt->setText(nui);
}
