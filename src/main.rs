
#[macro_use] extern crate custom_derive;
#[macro_use] extern crate newtype_derive;
use std::rc::Rc;
use std::sync::Arc;

use rug::{Float,Assign};
use kiss3d::camera::FirstPerson;
use kiss3d::nalgebra::{Vector3, UnitQuaternion, Point3};
use kiss3d::scene::SceneNode;
use kiss3d::window::Window;
use kiss3d::light::Light;

use stopwatch::Stopwatch;
pub mod simulation_float;
pub mod vec3;
pub mod pointmass;
pub mod simulation;
pub mod rope;
pub use vec3::*;
pub use simulation_float::*;
pub use pointmass::*;
pub use simulation::*;
pub use rope::*;

fn main() {
    let _config = SimulationConfig{
        g: SimulationFloat::new(6.67408e-11_f64),
        do_visuals: true,
    };
    //run_earth_moon_graphic();
    // run_simple_rotation();
    run_simple_rope_graphic();
}

fn run_simple_rope_graphic() {
    let mut window = Window::new("Trident Terminal Simulation");
    let mut camera = FirstPerson::new_with_frustrum(
        1.0, // fov
        0.00001, // znear
        1e+3_f32, // zfar
        Point3::new(0.0,6378.140,0.05), // position 
        Point3::new(0e+3_f32,6378.140_f32,0.0)); // pointing at

    let mut sim = Simulation::default();
    let mut start_model = window.add_cube(0.001, 0.001, 0.001);
    start_model.set_color(1.0,0.0,0.0);
    let mut end_model = window.add_cube(0.001, 0.001, 0.001);
    end_model.set_color(1.0,0.0,1.0);

    let start = PointMass::new(
            Vec3::new(10_f64, 6378140_f64 ,0_f64), // pos (at earth surface)
            Vec3::new(0_f32,0_f32,0_f32), // vel
            SimulationFloat::new(1_f64), // mass
            Some(start_model), // model
        );
    let end = PointMass::new(
            Vec3::new(-10_f64, 6378140_f64,0_f64), // pos (at earth surface)
            Vec3::new(0_f32,0_f32,0_f32), // vel
            SimulationFloat::new(1_f64), // mass
            Some(end_model), // model
        );
    sim.point_masses.push(start);
    sim.point_masses.push(end);
    sim.point_masses.push(PointMass::new_earth(Some(&mut window)));
    //sim.point_masses.push(PointMass::new_moon(Some(&mut window)));

    //sim.point_masses[1].enable_position_recording();
    window.set_light(Light::StickToCamera);

    let simulation_jiffy_s = SimulationFloat::new(0.1_f64);
    let simulation_speed_multiplier = SimulationFloat::new(0.5_f64);
    let real_elapsed_time = Stopwatch::start_new();
    //let mut ft = Stopwatch::start_new();
    //let mut frame_time_box = Box::new(Stopwatch::start_new());
    ////let frame_time = Rc::new(ft);
    //let inner_frame_time = frame_time_box.clone();


    let criteria = move |simulation_elapsed_time_s : SimulationFloat| {
         // (inner_frame_time.elapsed_ms() < 1000) && 
        (simulation_elapsed_time_s *
            SimulationFloat::new(1000_f32)) < simulation_speed_multiplier.clone() * SimulationFloat::new(real_elapsed_time.elapsed_ms())
        };

    while window.render_with_camera(&mut camera) {
        // frame_time_box.as_mut().restart();
        // dbg!(frame_time_box.elapsed_ms());
        sim.run_till(&simulation_jiffy_s, criteria.clone());
        for point in &mut sim.point_masses {
            point.update_scene_node_position();
        }

        sim.draw_historical_positions(&mut window);
        

    }


}

fn run_earth_moon_graphic() {
    let mut window = Window::new("Trident Terminal Simulation");
    let mut camera = FirstPerson::new_with_frustrum(
        1.0, // fov
        1.0, // znear
        1e+9_f32, // zfar
        Point3::new(0.0,0.0,700e+3_f32), // position 
        Point3::new(0e+3_f32,0.0,0.0)); // pointing at

    let mut sim = Simulation::default();
    sim.point_masses.push(PointMass::new_earth(Some(&mut window)));
    sim.point_masses.push(PointMass::new_moon(Some(&mut window)));
    //sim.point_masses.push(PointMass::new_moon(Some(&mut window)));

    sim.point_masses[1].enable_position_recording();
    window.set_light(Light::StickToCamera);

    let simulation_jiffy_s = SimulationFloat::new(50_f32);
    let simulation_speed_multiplier = SimulationFloat::new(4248_f32 * 3.0_f32);
    let real_elapsed_time = Stopwatch::start_new();
    //let mut ft = Stopwatch::start_new();
    //let mut frame_time_box = Box::new(Stopwatch::start_new());
    ////let frame_time = Rc::new(ft);
    //let inner_frame_time = frame_time_box.clone();


    let criteria = move |simulation_elapsed_time_s : SimulationFloat| {
         // (inner_frame_time.elapsed_ms() < 1000) && 
        (simulation_elapsed_time_s *
            SimulationFloat::new(1000_f32)) < simulation_speed_multiplier.clone() * SimulationFloat::new(real_elapsed_time.elapsed_ms())
        };

    while window.render_with_camera(&mut camera) {
        // frame_time_box.as_mut().restart();
        // dbg!(frame_time_box.elapsed_ms());
        sim.run_till(&simulation_jiffy_s, criteria.clone());
        for point in &mut sim.point_masses {
            point.update_scene_node_position();
        }

        sim.draw_historical_positions(&mut window);
        

    }
}
fn run_simple_rotation() {
    let mut window = Window::new("Trident Terminal Simulation");
    let mut camera = FirstPerson::new_with_frustrum(
        1.0, // fov
        1.0, // znear
        1e+9_f32, // zfar
        Point3::new(0.0,0.0,20_f32), // position 
        Point3::new(0e+3_f32,0.0,0.0)); // pointing at

    let mut sim = Simulation::default();
    let mut model1 = window.add_sphere(0.1_f32);
    let mut model2 = window.add_sphere(0.1_f32);
    model1.set_color(1.0,0.0,0.0);
    model2.set_color(1.0,0.0,1.0);
    let mass= 1e13_f64;
    sim.point_masses.push(PointMass::new_gravity(
            Vec3::new(-3000_f32,0_f32,30_f32), // pos (at earth surface)
            Vec3::new(0_f32,0.1_f32,0_f32), // vel
            SimulationFloat::new(mass), // mass
            Some(model1), // model
        ));
    sim.point_masses.push(PointMass::new_gravity(
            Vec3::new(3000_f32,0_f32,0_f32), // pos (at earth surface)
            Vec3::new(0_f32,-0.1_f32,0_f32), // vel
            SimulationFloat::new(mass), // mass
            Some(model2), // model
        ));
    //sim.point_masses.push(PointMass::new_moon(Some(&mut window)));

    window.set_light(Light::StickToCamera);

    let simulation_jiffy_s = SimulationFloat::new(5_f32);
    let simulation_speed_multiplier = SimulationFloat::new(4248_f32 * 3.0_f32);
    let real_elapsed_time = Stopwatch::start_new();

    let criteria = move |simulation_elapsed_time_s : SimulationFloat| {(
            simulation_elapsed_time_s *
            SimulationFloat::new(1000_f32)) < simulation_speed_multiplier.clone() * SimulationFloat::new(real_elapsed_time.elapsed_ms())
        };

    while window.render_with_camera(&mut camera) {
        sim.run_till(&simulation_jiffy_s, criteria.clone());
        for (index, point) in &mut sim.point_masses.iter_mut().enumerate() {
            point.update_scene_node_position();
        }

        

    }
}


#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn big_float_vector2() {
        let a = Vec3::new(-5.0_f32,1.0_f32,2.0_f32);
        dbg!(a);
        
    }
}
