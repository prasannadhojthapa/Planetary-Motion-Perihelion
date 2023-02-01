/* Physics- 215 Computational Math/Physics
//Lab 6
// Planetary Motion Perihelion 
*/

void Planetary_Motion_Perihelion(float pi = TMath::Pi())
{
	
	const int n = 12000;
	float dt = 0.01;
	//float v0 = 5; 
	float x0 = 1.0; 
	float t[n], x[n], y[n], v_x[n],v_y[n], r_1[n], r_2[n],r[n],v[n];
	float kinetic_by_mass[n], potential_by_mass[n] , energy_total[n];
	float angular_momentum_total_by_mass[n], d_A[n], theta_1[n], theta[n];
	float approximation_value = 0.1;
	float a_value, b_value, a_value_temp, b_value_temp;
	int b_value_coordinate;
	
	//coordinate of the sun
	a_value_temp = b_value_temp =  0.;

	
	//Intial values declaration starts here:
	//This is the extrensity
	float E = 0.206;
	a_value = 0.39; // AU unit
	b_value = a_value*sqrt(1 - E*E);
	r_1[0] = a_value*(1+E);
	r_2[0] = sqrt(E*a_value*E*a_value + b_value*b_value);
	v_x[0] = 0;
	//For the mercury
	v_y[0] =  sqrt( (4*pi*pi*(1-E)) / (a_value*(1+E)) ); 

	x[0] = x0;
	y[0] = 0;
	
	v[0] = sqrt(v_x[0]*v_x[0] +v_y[0]*v_y[0] );
	r[0] = sqrt(x[0]*x[0] + y[0]*y[0]);
	
	//Initializing the first index of Kinetic, Potential and Total energy array.
	kinetic_by_mass[0] = 0.5 * v[0] * v[0];
	potential_by_mass[0] = ( 4*pi*pi) / r[0]; 
	energy_total[0] = kinetic_by_mass[0] + potential_by_mass[0]; 
	
	//for the angle
	theta_1[0]=0.;
	theta[0] = acos(x[0]/r[0]);
	d_A[0] = 0;
	
	angular_momentum_total_by_mass[0] = x[0]*v_y[0] - y[0]*v_x[0];
	
	cnvs_planetary_motion = new TCanvas("cnvs_planetary_motion", "cnvs_planetary_motion",10,10,1300, 800);
	cnvs_planetary_motion -> Draw();

	for(int i=0; i<n-1; i++)
	{
		
		t[i+1] = t[i] + dt;
		
		v_x[i+1] = v_x[i] - (4*pi*pi*x[i]*dt)/(r[i]*r[i]*r[i]);
		v_y[i+1] = v_y[i] - (4*pi*pi*y[i]*dt)/(r[i]*r[i]*r[i]);

		x[i+1] = x[i] + v_x[i+1]*dt;
		y[i+1] = y[i] + v_y[i+1]*dt;
		
		r[i+1] = sqrt(x[i+1]*x[i+1] + y[i+1]*y[i+1]); 
		v[i+1] = sqrt(v_x[i+1]*v_x[i+1] +v_y[i+1]*v_y[i+1] );
		
		//For the energ y conservation
		kinetic_by_mass[i+1] = 0.5 * v[i+1] * v[i+1];
		potential_by_mass[i+1] = -(4*pi*pi) / r[i+1]; 
		energy_total[i+1] = kinetic_by_mass[i+1] + potential_by_mass[i+1]; 
		
		//For the angular momentum conservation
		angular_momentum_total_by_mass[i+1] = x[i+1]*v_y[i+1] - y[i+1]*v_x[i+1];
		
		//area of triangle:
		//for the angle
		theta[i+1] = acos(x[i+1]/r[i+1]);
//		d_A[i+1] = 0.5*r[i+1]*r[i+1]*fabs(theta[i+1]-theta[i]);
		d_A[i+1] = 0.5*r[i+1]*sqrt((x[i+1]-x[i])*(x[i+1]-x[i]) + (y[i+1]-y[i])*(y[i+1]-y[i]));

		
		
		if( a_value < fabs(x[i+1]))
		{
			a_value = fabs(x[i+1]); 
			a_value_temp =  fabs(x[i+1]) - fabs(x[b_value_coordinate]);

		}
		if( b_value < fabs(y[i+1]) )
		{
			b_value = fabs(y[i+1]); 
			b_value_temp = fabs(y[i+1]);
			b_value_coordinate = i+1;			
			
		}
	
	
	}
	
	gr_planetary_motion_through_time = new TGraph(n,x,y);
	gr_planetary_motion_through_time -> SetTitle("Planetary Motion");
	gr_planetary_motion_through_time -> Draw("AL");
	
	
	circle_sun = new TEllipse(0,0,0.3,0.3);
	circle_sun->SetFillColor(5);
	circle_sun->SetFillStyle(1001);
	circle_sun-> SetLineColor(5);
	circle_sun->Draw("same");
	
	circle_planet1 = new TEllipse(x[n-1],y[n-1],0.1,0.1);
	circle_planet1->SetFillColor(1);
	circle_planet1->SetFillStyle(3008);
	circle_planet1->Draw("same");
	
	//Line for the a value; the semi-major axis
	line_for_a_value = new TLine(x[b_value_coordinate],0,-a_value, 0);
	line_for_a_value -> SetLineColor(2);
	line_for_a_value -> SetLineStyle(5);
	line_for_a_value -> Draw("same");
	
	//Line for the b value; the semi-minor axis
	
	line_for_b_value = new TLine(x[b_value_coordinate],0,x[b_value_coordinate],b_value_temp);
	line_for_b_value -> SetLineColor(4);
	line_for_b_value -> SetLineStyle(5);
	line_for_b_value -> Draw("same");
	
	char name1[100];
	char name2[100];
	
	
	
	auto leg_planetary_motion = new TLegend(0.7,0.85,0.9,1);
	sprintf(name1,"Semi-major axis: %f",a_value_temp);
	sprintf(name2,"Semi-minor axis: %f",b_value_temp);
	
	
	//Dummy line for the legend section [Planet orbit]
	line_for_the_orbit = new TLine();
	line_for_the_orbit-> SetLineColor(1);
	line_for_the_orbit-> Draw("same");
	
	leg_planetary_motion -> AddEntry(line_for_a_value,name1);
	leg_planetary_motion -> AddEntry(line_for_b_value,name2);
	leg_planetary_motion -> AddEntry(line_for_the_orbit,"Planet orbit");
	leg_planetary_motion -> AddEntry(circle_sun,"Sun");
	leg_planetary_motion -> AddEntry(circle_planet1,"Planet");
	leg_planetary_motion -> SetBorderSize(1);
	leg_planetary_motion -> Draw("same");
		
	//Save file
	cnvs_planetary_motion -> SaveAs("C:/root_v5.34.36/plot_planetary_motion_orbit_graph_Eliptical_position_1AU_initialvelocity_5AU.gif");

	
	cnvs_energy_Energy = new TCanvas("cnvs_planetary_motion_energy","cnvs_planetary_motion_energy",10,10,1400,800);
	cnvs_energy_Energy->Draw();
	cnvs_energy_Energy -> Divide(3,1);
	cnvs_energy_Energy->cd(1);
	
	
	gr_energy_kinetic = new TGraph(n,t, kinetic_by_mass);
	gr_energy_kinetic->SetTitle("Kinetic Energy vs Time");
	gr_energy_kinetic->GetYaxis()->SetTitle("Kinetic");
	gr_energy_kinetic->GetXaxis()->SetTitle("Time, seconds");
	//This sets our desired y-axis range
	gr_energy_kinetic->GetYaxis()->SetRangeUser(0,50);
	gr_energy_kinetic-> SetLineStyle(20);
	gr_energy_kinetic-> SetLineColor(1);
	gr_energy_kinetic-> Draw("AL");
	
	
	cnvs_energy_Energy->cd(2);
	
	gr_energy_potential = new TGraph(n,t, potential_by_mass);
	gr_energy_potential->SetTitle("Potential Energy vs Time");
	gr_energy_potential->GetYaxis()->SetTitle("Potential");
	//This sets our desired y-axis range
	gr_energy_potential->GetYaxis()->SetRangeUser(-50,0);
	//continued code
	gr_energy_potential->GetXaxis()->SetTitle("Time, seconds");
	gr_energy_potential-> SetLineStyle(20);
	gr_energy_potential-> SetLineColor(1);
	gr_energy_potential-> Draw("AL");
	
	cnvs_energy_Energy->cd(3);
	
	gr_energy_total = new TGraph(n,t, energy_total);
	gr_energy_total->SetTitle("Total Energy vs Time");
	gr_energy_total->GetYaxis()->SetTitle("Total Energy");
	gr_energy_total->GetXaxis()->SetTitle("Time, seconds");
	//This sets our desired y-axis range
	//gr_energy_total->GetYaxis()->SetRangeUser(0,20);
	//continued code
	gr_energy_total-> SetLineStyle(20);
	gr_energy_total-> SetLineColor(1);
	gr_energy_total-> Draw("AP");
	
	//Save file
	cnvs_energy_Energy -> SaveAs("C:/root_v5.34.36/plot_planetary_motion_all_energy_graph_Eliptical_position_1AU_initialvelocity_5AU.gif");

	
	cnvs_planetary_angular_momentum_graph_and_area = new TCanvas("cnvs_planetary_angular_momentum_graph","cnvs_planetary_angular_momentum_graph",10,10,1400,800);
	cnvs_planetary_angular_momentum_graph_and_area->Draw();
	cnvs_planetary_angular_momentum_graph_and_area -> Divide(2,1);
	
	cnvs_planetary_angular_momentum_graph_and_area->cd(1);
	
	gr_angular_momentum = new TGraph(n,t, angular_momentum_total_by_mass);
	gr_angular_momentum->SetTitle("Angular Momentum vs Time");
	gr_angular_momentum->GetYaxis()->SetTitle("Angular momentum");
	gr_angular_momentum->GetXaxis()->SetTitle("Time, seconds");
	//This sets our desired y-axis range
	gr_angular_momentum->GetYaxis()->SetRangeUser(0,20);
	//continued code
	gr_angular_momentum-> SetLineStyle(20);
	gr_angular_momentum-> SetLineColor(1);
	gr_angular_momentum-> Draw("AL");
	
	cnvs_planetary_angular_momentum_graph_and_area->cd(2);
	gr_d_A = new TGraph(n,t,d_A);
	gr_d_A->SetTitle("dA vs Time");
	gr_d_A->GetYaxis()->SetTitle("dA");
	gr_d_A->GetXaxis()->SetTitle("Time, seconds");
	gr_d_A-> SetLineStyle(20);
	gr_d_A-> SetLineColor(1);
	gr_d_A-> Draw("AL");
	
	cnvs_planetary_angular_momentum_graph_and_area -> SaveAs("C:/root_v5.34.36/plot_planetary_motion_angular_momentum_and_area_graph_Eliptical_position_1AU_initialvelocity_5AU.gif");
	
}
