// Copyrigth 2008-2017 by Edmanuel Torres A. (eetorres@gmail.com)
//
// This software is  free  software;  you  can  redistribute it and/or
// modify it  under  the  terms  of  the   GNU  Library General Public
// License  as  published  by  the  Free  Software Foundation;  either
// version 2 of the License,  or  (at your option)  any later version.
// or even much better the FLTK license wich allow you static linking.
//
// This  library  is  distributed  in the hope that it will be useful,
// but  WITHOUT ANY WARRANTY;  without  even  the  implied warranty of
// MERCHANTABILITY  or FITNESS FOR A PARTICULAR PURPOSE.   See the GNU
// Library General Public License for more details.
//
// You should have  received a copy  of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
// The above copyright  notice  and  this  permission notice should be
// included in all  the copies  or  modifications on this code. I will
// be very grateful  to  receive your feedbacks.   For any suggestion,
// modification or bugs. Don't hesitate to contact me. Get the lastest
// version at https://github.com/eetorres/mdargon

// C++ libraries
#include <unistd.h>

// OpenGL
#include <GL/gl.h>
#include <GL/glut.h>

// MD header
#include <md3d.h>

int sample_energy = 10;
int sample_posicion = 100;
int steps = 0;

void write_coordinates(void){
    // file to store the coordinates
    argon_xyz.open("argon_ini.xyz");
    for (int i = 0; i < N; i++){
        // save the coordiantes in the file
        argon_xyz<<r[i][0]<<"  ";
        argon_xyz<<r[i][1]<<"  ";
        argon_xyz<<r[i][2]<<"  "<<endl;
    }
    argon_xyz.close();
}


// Compute the trajectories for 100000 time steps
void motion(void){
    //
    velocityVerlet(dt);
    //cout <<" T: "<< instantTemperature();
    //cout <<" U: "<< PotentialEnergy();
    //cout <<" K: "<< KineticEnergy() << '\n';
    //cout <<" Total energy : "<< KineticEnergy()+PotentialEnergy() << " (per particle)\n";
    steps++;
    if((steps%sample_energy)==0){
        argon_energy<<steps<<"  "<<potentialEnergy()<<"  "<<kineticEnergy()<<"  "<<kineticEnergy()+potentialEnergy()<<endl;
        cout<<steps<<"  "<<potentialEnergy()<<"  "<<kineticEnergy()<<"  "<<kineticEnergy()+potentialEnergy()<<endl;
    }
}

// Display the simulation
void point_draw(void) {
    glPushMatrix();
    // show the cubic box
    glColor3f(0, 1, 0);
    glutWireCube((GLfloat)L);
    // Show the position of the  particles
    glPointSize(1);
    glColor3f(1, 0, 0);
    glBegin(GL_POINTS);
    for(int i=0; i<N; i++)
       glVertex3f((GLfloat)r[i][0], (GLfloat)r[i][1], (GLfloat)r[i][2]);
    glEnd();
    glPopMatrix();
    glPopMatrix();
}

// OpenGL configuration
void display(void) {
    float lado = 1.5*L_2;
    glClearColor (0.0,0.0,0.0,1.0);
    glLoadIdentity();
    glOrtho(lado,-lado,lado,-lado,lado,-lado);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glRotatef(20,1,0,0);
    glRotatef(10,0,1,0);
    //
    point_draw();
    motion();
    //
    glutSwapBuffers();
    // delay of 50 microseconds
    usleep(500);
}

// Resize windows
void reshape (int w, int h) {
    glViewport (0, 0, (GLsizei)w, (GLsizei)h);
    glLoadIdentity ();
    glMatrixMode (GL_MODELVIEW);
}

/// Main function
int main (int argc, char **argv) {
    glutInit (&argc, argv);
    // initial condisitons
    initialize();
    // file to store the energies
    argon_energy.open("argon_energy.dat");
    //
    cout<<"Box side \t= \t"<<L<<endl;
    cout<<"Volume \t\t= \t"<<L*L*L<<endl;
    cout<<"Density \t= \t"<<rho<<endl;
    cout<<"Lattice \t= \t"<<lattice_contant<<endl<<endl;
    cout<<"N      U        K        E"<<endl;
    //
    // mode of the display
    glutInitDisplayMode (GLUT_DOUBLE); //set up the double buffering
    // dimensiones de la venana
    glutInitWindowSize (500, 500);
    // position of the window
    glutInitWindowPosition (100, 100);
    // create a OpenGL window instance
    glutCreateWindow ("Molecular dynamics of argon");
    // glut functions
    glutDisplayFunc (display);
    glutIdleFunc (display);
    glutReshapeFunc (reshape);
    glutMainLoop ();
    // end of the MD simulation
    return 0;
}

//
