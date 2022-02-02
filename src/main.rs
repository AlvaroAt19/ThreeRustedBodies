extern crate ndarray;
use ndarray::prelude::*;
use poloto::prelude::*;


fn main() {
    let mut values: Array<f64,Ix1> = arr1(&[-1.0,0.322184765624991,0.0,0.647989160156249,1.0,0.322184765624991,0.0,0.647989160156249,0.0,-2.0*0.322184765624991,0.0,-2.0*0.647989160156249]);
    
    let mut rkx_pnts1:Vec<f64> = [-1.0].to_vec();
    let mut rky_pnts1:Vec<f64> = [0.0].to_vec();


    let mut rkx_pnts2:Vec<f64> = [1.0].to_vec();
    let mut rky_pnts2:Vec<f64> = [0.0].to_vec();

    let mut rkx_pnts3:Vec<f64> = [0.0].to_vec();
    let mut rky_pnts3:Vec<f64> = [0.0].to_vec();


    let h:f64 = 0.0001;

    let mut t:f64 = 0.0;

    while t < 51.3958{

        let k1 = u(&values,t);

        let vk2 = &values + ((h*&k1)/2.0); //2° step array

        let k2 = u(&vk2, (t+h)/2.0);
        
        let vk3 = &values + ((h*&k2)/2.0); //3° step array

        let k3 = u(&vk3, (t+h)/2.0);

        let vk4 = &values + (h*&k3); //4° step array

        let k4 = u(&vk4, t+h);

        values = &values + h*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

        let copy = values.clone();

        t = t + h;

        rkx_pnts1.push(copy[0]);
        rky_pnts1.push(copy[2]);
        rkx_pnts2.push(copy[4]);
        rky_pnts2.push(copy[6]);
        rkx_pnts3.push(copy[8]);
        rky_pnts3.push(copy[10]);
    };

    let mut plotter = poloto::plot(
        "Trajectory of 1°st body",
        formatm!("This is the {} label", 'x'),
        "This is the y label",
        f64::default_ctx(),
        f64::default_ctx(),
    );

    let mut coord = Vec::new();

    for i in 0..(5130){
        coord.push([rkx_pnts1[i*100],rky_pnts1[i*100]]);
    }

    plotter.line("tragectory", coord);

    println!("{}", poloto::disp(|a| plotter.simple_theme(a)));
}


fn u(r:&Array<f64,Ix1>, _t:f64)->Array<f64,Ix1>{
    let g = 1.0;  //Gravity
    let m1 = 1.0; //Mass of 1° body
    let m2 = 1.0; //Mass of 2° body
    let m3 = 1.0; //Mass of 3° body

    
    let x1 = r[0];  //X Position of 1° body
    let v_x1 = r[1];//X Velocity of 1° body
    let y1 = r[2];  //Y Position of 1° body
    let v_y1 = r[3];//Y Velocity of 1° body

    let x2 = r[4];       
    let v_x2 = r[5];
    let y2 = r[6];
    let v_y2 = r[7];
    
    let x3 = r[8];       
    let v_x3 = r[9];
    let y3 = r[10];
    let v_y3 = r[11];

    //x and y aceleration

    let dvx1:f64 = ((-(g*m2)*(x1 - x2))/(((x1-x2).powf(2.0) + (y1-y2).powf(2.0)).powf(1.5))) + ((-(g*m3)*(x1 - x3))/(((x1-x3).powf(2.0) + (y1-y3).powf(2.0)).powf(1.5))); 
    let dvy1:f64 = ((-(g*m2)*(y1 - y2))/(((x1-x2).powf(2.0) + (y1-y2).powf(2.0)).powf(1.5))) + ((-(g*m3)*(y1 - y3))/(((x1-x3).powf(2.0) + (y1-y3).powf(2.0)).powf(1.5)));    
    let dvx2:f64 = ((-(g*m1)*(x2 - x1))/(((x1-x2).powf(2.0) + (y1-y2).powf(2.0)).powf(1.5))) + ((-(g*m3)*(x2 - x3))/(((x2-x3).powf(2.0) + (y2-y3).powf(2.0)).powf(1.5)));  
    let dvy2:f64 = ((-(g*m1)*(y2 - y1))/(((x1-x2).powf(2.0) + (y1-y2).powf(2.0)).powf(1.5))) + ((-(g*m3)*(y2 - y3))/(((x2-x3).powf(2.0) + (y2-y3).powf(2.0)).powf(1.5)));
    let dvx3:f64 = ((-(g*m1)*(x3 - x1))/(((x1-x3).powf(2.0) + (y1-y3).powf(2.0)).powf(1.5))) + ((-(g*m2)*(x3 - x2))/(((x2-x3).powf(2.0) + (y2-y3).powf(2.0)).powf(1.5)));  
    let dvy3:f64 = ((-(g*m1)*(y3 - y1))/(((x1-x3).powf(2.0) + (y1-y3).powf(2.0)).powf(1.5))) + ((-(g*m2)*(y3 - y2))/(((x2-x3).powf(2.0) + (y2-y3).powf(2.0)).powf(1.5))); 


    let resp:Array<f64,Ix1> = arr1(&[v_x1, dvx1, v_y1, dvy1, v_x2, dvx2, v_y2, dvy2,v_x3,dvx3,v_y3,dvy3]);

    return resp;
}