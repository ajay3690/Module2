#include <Arduino.h>
#ifdef ESP32
  #include <WiFi.h>
  #include <AsyncTCP.h>
#else
  #include <ESP8266WiFi.h>
  #include <ESPAsyncTCP.h>:
#endif
#include <ESPAsyncWebServer.h>
#include"matfun.h"

AsyncWebServer server(80);

const char* ssid = "Ajay";
const char* password = "Ajay2610";

const char* input_parameter00 = "input00";
const char* input_parameter01 = "input01";
const char* input_parameter02 = "input02";
const char* input_parameter10 = "input10";
const char* input_parameter11 = "input11";
const char* input_parameter12 = "input12";
const char* input_parameter20 = "input20";
const char* input_parameter21 = "input21";
const char* input_parameter22 = "input22";
const char* matrix1[3]={input_parameter00,input_parameter01,input_parameter02};     // matrix for vertex B
const char* matrix2[3]={input_parameter10,input_parameter11,input_parameter12};     // matrix for vertex C
const char* matrix3[3]={input_parameter20,input_parameter21,input_parameter22};     // matrix for vertex D

const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML><html><head>
    <title>TRIANGLE</title>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
      html {font-family: Times New Roman; display: inline-block; text-align: center;}
      h2 {font-size: 2.0rem; color: blue;}
    </style> 
    </head><body>
    <h2>Collinearity of points using gram-schmit</h2> 
    <p>Enter the values of points v1, v2 and v3
    <form action="/get">
      Enter the values of Point v1: <input type="number" name="input00"> <input type="number" name="input01"><input type="number" name="input02"><br><br>
      Enter the values of Point v2: <input type="number" name="input10"> <input type="number" name="input11"><input type="number" name="input12"><br><br>
      Enter the values of Point v3: <input type="number" name="input20"> <input type="number" name="input21"><input type="number" name="input22"><br><br> 

      <input type="submit" value="Submit">
    </form><br>
  </body></html>)rawliteral";

void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void setup() {
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  if (WiFi.waitForConnectResult() != WL_CONNECTED) {
    Serial.println("Connecting...");
    return;
  }
  Serial.println();
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

server.on("/get", HTTP_GET, [] (AsyncWebServerRequest *request) {


double **v1,**v2,**v3,**A,**B,**C,**D,**M,**F,**G; //declaring matrices names


v1=load_ser(request,matrix1,3);
v2=load_ser(request,matrix2,3);
v3=load_ser(request,matrix3,3);

/*printf("\n v1 \n");
printMat(v1,3,1);
printf("\n v2 \n");
printMat(v2,3,1);
printf("\n v3 \n");
printMat(v3,3,1);*/

double vn1,vn2,vn3,**b1,**b11;

vn1 = norm(v1);
//printf("\n norm of v1 vn1 = %lf\n\n ",vn1);
A=Matscale(v1,3,1,1/vn1); //A= v1/vn1
//printf("\n A =\n");
//printMat(A,3,1);

// B

b1=tr_mul_mat(v2,A,A);// (v2.T @A )* A
b11=Matsub(v2,b1,3,1); // v2- b111
vn2=norm(b11); // norm b11
//printf("\n norm of b1 = %lf \n",vn2);

B= Matscale(b11,3,1,1/vn2);
//printf("\n B\n");
//printMat(B,3,1);

// C
double **v3_a,**v3_b,**c1,**c11;
v3_a=tr_mul_mat(v3,A,A);//(v3.T @ A )* A
v3_b=tr_mul_mat(v3,B,B);//(v3.T @  B)* B

c11=Matsub(v3,v3_a,3,1);
c1=Matsub(c11,v3_b,3,1);
vn3=norm(c1);
//printf("\n norm of c1 vn3 = %lf\n\n ",vn3);
C= Matscale(c1,3,1,1/vn3);

//printf("\n C \n");
//printMat(C,3,1);

double **ABC, **ABC_T,**I;
ABC=p_m(A,B,C,3,3); // [A B C]
//printf("\n [ A B C ] \n");
//printMat(ABC,3,3);

ABC_T=transposeMat(ABC,3,3);
//printf("\n [ A B C ].T \n");
//printMat(ABC_T,3,3);

I= Matmul(ABC,ABC_T,3,3,3); // ABC*ABC.T
//printf("\n Identity matrix \n");
//printMat(I,3,3);

double u1,u2,u3;
u1=norm(Matcol(I,3,0));// norm of 0th column of ABC
u2=norm(Matcol(I,3,1));// norm of 1th column of ABC
u3=norm(Matcol(I,3,2));// norm of 2nd column of ABC
//printf("\nu1= %lf, u2= %lf, u3=%lf\n",u1,u2,u3);

String response;

if(round(u1)==round(u2)==round(u3)){
	D=Matadd(A,Matadd(B,C,3,1),3,1); // D= A+B+C
	//printf("\n D\n");
	//printMat(D,3,1);
        double D_norm,C_norm,B_norm,A_norm;
	D_norm=norm(D);
	C_norm=norm(C);
	B_norm=norm(B);
	A_norm=norm(A);
	double theta1,theta2,theta3;
	theta1=Angle(A,D,A_norm,D_norm);
	//printf("\nTheta1=%lf\n",theta1);
	theta2=Angle(B,D,B_norm,D_norm);
	//printf("\nTheta2=%lf\n",theta2);
	theta3=Angle(C,D,C_norm,D_norm);
	//printf("\nTheta3=%lf\n",theta3);
	response +="<p> A+B+c equally inclined to A,B,C</p>";
	//printf("\n A+B+c equally inclined to A,B,C\n");
}
else{
	//printf("\n A+B+C are not equally inclined to A, B, C\n");
	response +="<p> A+B+c are not equally inclined to A,B,C</p>";
}


	response += "<br><a href=\"/\">Return to Home Page</a>";
    // Send the HTML response with dynamic content
    request->send(200, "text/html", response);
});
  server.onNotFound(notFound);
  server.begin();
}
void loop() { 
}
