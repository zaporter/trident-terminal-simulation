use kiss3d::{scene::SceneNode, nalgebra::Vector3};
use kiss3d::window::Window;

use std::fmt::Debug;
use std::fmt;
use uuid::Uuid;
use std::collections::HashMap;
use crate::*;

#[derive( Clone)]
pub struct PointMass {
    pub pos: Vec3,
    pub vel : Vec3,
    pub mass: SimulationFloat,
    pub scene_node : Option<SceneNode>,
    pub applies_gravity:bool,
    pub temp_K : SimulationFloat,
    external_forces : HashMap<Uuid, Vec3>,
    past_positions : Vec<Vec3>,
    should_record_positions : bool,
}
impl PointMass {
    pub fn new(pos :Vec3, vel:Vec3, mass:SimulationFloat, scene_node : Option<SceneNode>) ->PointMass {
        PointMass { 
            pos, 
            vel, 
            mass, 
            scene_node,
            applies_gravity:false,
            temp_K : SimulationFloat::default(),
            past_positions : Vec::new(),
            should_record_positions : false,
            external_forces : HashMap::new(),
        }
    }
    pub fn new_gravity(pos :Vec3, vel:Vec3, mass:SimulationFloat, scene_node : Option<SceneNode>) ->PointMass {
        PointMass { 
            pos, 
            vel, 
            mass, 
            scene_node,
            applies_gravity:true,
            temp_K : SimulationFloat::default(),
            past_positions : Vec::new(),
            should_record_positions : false,
            external_forces : HashMap::new(),
        }
    }
    pub fn new_earth(window:Option<&mut Window>) -> PointMass{
        // https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
        //let earth_radius_km = 3671_f32;
        let earth_radius_km = -6378.137_f32;
        let earth_mass_kg = SimulationFloat::new(5.9724e+24_f64);
        let mut earth_model : Option<SceneNode> = None;
        if let Some(window) = window {
            let mut earth_model_unwrapped = window.add_sphere(earth_radius_km);
            earth_model_unwrapped.set_color(0.0,0.4,0.7);
            earth_model = Some(earth_model_unwrapped);
        }
        let earth_point_mass = PointMass::new_gravity(
                Vec3::new(0_f32,0_f32,0_f32), // pos
                Vec3::new(0_f32,0_f32,0_f32), // vel
                earth_mass_kg, // mass
                earth_model, // model
            );
        earth_point_mass

    }
    pub fn new_moon(window: Option<&mut Window>) -> PointMass {
        // https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
        let moon_radius_km = 1738_f32;
        let moon_mass_kg = SimulationFloat::new(0.07346e+24_f64);
        let mut moon_model : Option<SceneNode> = None;
        if let Some(window) = window {
            let mut moon_model_unwrapped = window.add_sphere(moon_radius_km);
            moon_model_unwrapped.set_color(0.9,0.9,0.9);
            moon_model = Some(moon_model_unwrapped);
        }
        let mut moon_point_mass = PointMass::new_gravity(
                Vec3::new(333344663_f64,0_f64,145458387_f64), // pos
                Vec3::new(0_f32,1082_f32,0_f32), // vel
                moon_mass_kg, // mass
                moon_model, // model
            );
        moon_point_mass
    }
    
    pub fn calculate_gravitational_force(&self, gravity_applying_masses : &Vec<&PointMass>) -> Vec3 {
        let mut force = Vec3::new(0_f32, 0_f32, 0_f32);
        let g_m3_kg_s = SimulationFloat::new(6.673889e-11_f64);
        for gravity_point_mass in gravity_applying_masses {
            let difference = gravity_point_mass.pos.clone() - self.pos.clone();

            let mut distance_squared =  difference.modulus();
            distance_squared.square_mut();
            if distance_squared == SimulationFloat::new(0_f32) {
                continue;
            }
            let difference_norm = difference.unit();
            let product_of_masses = gravity_point_mass.mass.clone() * self.mass.clone();

            let grav_force_magnitude = g_m3_kg_s.clone() * (product_of_masses / distance_squared);
            let grav_force = difference_norm * grav_force_magnitude;

            force += grav_force;
        }
        force
    }
    pub fn print_total_energy_j(&self, gravity_applying_masses : Vec<&PointMass>) {
        dbg!(self.calculate_kinetic_energy_j());
            dbg!(self.calculate_potential_spring_energy_j());
            dbg!(self.calculate_gravitational_potential_energy_j(gravity_applying_masses)); 
            dbg!(self.calculate_heat_energy_j());
    }
    pub fn enable_position_recording(&mut self) {
        self.should_record_positions = true;
    }
    pub fn disable_position_recording(&mut self) {
        self.should_record_positions = false;
    }
        

    pub fn calculate_kinetic_energy_j(&self) -> SimulationFloat {
        let mut speed_squared = self.vel.modulus();
        speed_squared.square_mut();
        (self.mass.clone() * speed_squared ) / SimulationFloat::new(2_f32)
    }
    pub fn calculate_potential_spring_energy_j(&self) -> SimulationFloat {
        // -k (L-L0) / 2
        // divide this by two and do it on both sides. 
        SimulationFloat::default()
    }
    pub fn calculate_gravitational_potential_energy_j(&self, gravity_applying_masses : Vec<&PointMass>) -> SimulationFloat {
        let mut energy = SimulationFloat::default();
        let g_m3_kg_s = SimulationFloat::new(6.673889e-11_f64);
        for gravity_point_mass in gravity_applying_masses {
            let difference = gravity_point_mass.pos.clone() - self.pos.clone();

            let mut distance =  difference.modulus();
            if distance == SimulationFloat::new(0_f32) {
                continue;
            }
            let product_of_masses = gravity_point_mass.mass.clone() * self.mass.clone();

            let pot_energy_magnitude = g_m3_kg_s.clone() * (product_of_masses / distance);
            energy += pot_energy_magnitude;
        }
        energy * SimulationFloat::new(-1_f32)

    }
    pub fn calculate_se_energy_j(&self, gravity_applying_masses : Vec<&PointMass>) -> SimulationFloat {
        let mut energy = SimulationFloat::default();
        let g_m3_kg_s = SimulationFloat::new(6.673889e-11_f64);
        for gravity_point_mass in gravity_applying_masses {
            let difference = gravity_point_mass.pos.clone() - self.pos.clone();

            let mut distance =  difference.modulus();
            if distance == SimulationFloat::new(0_f32) {
                continue;
            }
            let product_of_masses = (gravity_point_mass.mass.clone() + self.mass.clone());

            let pot_energy_magnitude = g_m3_kg_s.clone() * (product_of_masses / distance);
            let mut es = (gravity_point_mass.vel.clone() - self.vel.clone()).modulus();
            es.square_mut();
            es = es/SimulationFloat::new(2.0);
            energy += (es - pot_energy_magnitude);
        }
        energy //* SimulationFloat::new(-1_f32)

    }
    pub fn calculate_heat_energy_j(&self) -> SimulationFloat {
        self.mass.clone() * self.temp_K.clone() // * SPECIFIC HEAT   
    }
    pub fn add_external_force(&mut self, force_id : &Uuid, force : Vec3) {
        self.external_forces.insert(*force_id, force);
    }
    pub fn remove_external_force(&mut self, force_id : &Uuid) -> Option<Vec3> {
        self.external_forces.remove(force_id)

    }
}

impl Simulatable for PointMass {
    fn tick(&mut self, dt_s : &SimulationFloat, gravity_applying_masses : &Vec<&PointMass>) {
        let mut force = Vec3::new(0_f32, 0_f32, 0_f32);
        force += self.calculate_gravitational_force(gravity_applying_masses);

        for (_, external_force) in &self.external_forces {
            force += external_force.clone();
        }

        let acceleration = force / self.mass.clone();

        self.vel += acceleration     * dt_s.clone();
        self.pos += self.vel.clone() * dt_s.clone();

        if self.should_record_positions {
            self.past_positions.push(self.pos.clone());
        }
    }
    fn calculate_momentum(&self) -> Vec3 {
        self.vel.clone() * self.mass.clone() 
    }
    fn calculate_total_energy_j(&self, gravity_applying_masses : Vec<&PointMass>) -> SimulationFloat {
        self.calculate_kinetic_energy_j() + 
            self.calculate_potential_spring_energy_j() +
            self.calculate_gravitational_potential_energy_j(gravity_applying_masses) + 
            self.calculate_heat_energy_j()
        //self.calculate_se_energy_j(gravity_applying_masses)
    }

    fn update_scene_node_position(&mut self){
        if let Some(node) = &mut self.scene_node{
            node.set_local_translation(kiss3d::nalgebra::Translation::from(Vector3::new(
                            self.pos.x.to_f32() / 1000.0,
                            self.pos.y.to_f32() / 1000.0,
                            self.pos.z.to_f32() / 1000.0,
                        )));
        }
    }
    fn draw_historical_positions(&self, window : &mut Window) {
        for line_pair in self.past_positions.windows(2) {
            window.draw_line(
                &Point3::new(
                    line_pair[0].x.to_f32()/1000.0,
                    line_pair[0].y.to_f32()/1000.0,
                    line_pair[0].z.to_f32()/1000.0,
                    ),
                &Point3::new(
                    line_pair[1].x.to_f32()/1000.0,
                    line_pair[1].y.to_f32()/1000.0,
                    line_pair[1].z.to_f32()/1000.0,
                    ),
                &Point3::new(
                    0.7,
                    0.3,
                    0.3,
                    )
                )
                
            }
    }

}

impl Debug for PointMass {
     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PointMass")
         .field("pos", &self.pos)
         .field("vel", &self.vel)
         .field("mass", &self.mass)
         .field("applies_gravity", &self.applies_gravity)
         .finish()
    }

}


#[cfg(test)]
mod tests{
    use crate::*;
    #[test]
    fn equal_grav_forces_sum_zero(){

        let planet = PointMass::new_gravity(
                Vec3::new(10e+6_f32,0_f32,0_f32), // pos
                Vec3::new(0_f32,0_f32,0_f32), // vel
                SimulationFloat::new(1e+9_f64), // mass
                None, // model
            );
        let moon = PointMass::new_gravity(
                Vec3::new(-10e+6_f32,0_f32,0_f32), // pos
                Vec3::new(0_f32,10000_f32,0_f32), // vel
                SimulationFloat::new(0.5e+9_f64), // mass
                None, // model
            );
        let p_clone = planet.clone();
        let m_clone = moon.clone();
        let grav_masses = vec![&p_clone, &m_clone];
        let planet_force = planet.calculate_gravitational_force(&grav_masses);
        let moon_force = moon.calculate_gravitational_force(&grav_masses);
        assert!((planet_force+moon_force).modulus() < SimulationFloat::new(10_f32));

    }
    #[test]
    fn real_grav_forces(){

        let planet = PointMass::new_gravity(
                Vec3::new(0_f32,0_f32,0_f32), // pos
                Vec3::new(0_f32,0_f32,0_f32), // vel
                SimulationFloat::new(5.972e+24_f64), // mass
                None, // model
            );
        let moon = PointMass::new_gravity(
                Vec3::new(0.3633e+9_f32,0_f32,0_f32), // pos
                Vec3::new(0_f32,10000_f32,0_f32), // vel
                SimulationFloat::new(0.07346e+24_f64), // mass
                None, // model
            );
        let p_clone = planet.clone();
        let m_clone = moon.clone();
        let grav_masses = vec![&p_clone, &m_clone];
        let planet_force = planet.calculate_gravitational_force(&grav_masses);
        let moon_force = moon.calculate_gravitational_force(&grav_masses);
        dbg!(planet_force);
        dbg!(moon_force);

    }

    #[test]
    fn calc_grav_earth_surface(){
        let earth = PointMass::new_earth(None);
        let onekilo = PointMass::new_gravity(
                Vec3::new(-6378137_f32,0_f32,0_f32), // pos (at earth surface)
                Vec3::new(0_f32,0_f32,0_f32), // vel
                SimulationFloat::new(1_f64), // mass
                None, // model
            );
        let force = onekilo.calculate_gravitational_force(&vec![&earth]).modulus();
        assert!(force > SimulationFloat::new(9.75_f32));
        assert!(force < SimulationFloat::new(9.85_f32));
    }
    #[test]
    fn calc_potential_energy_test(){
        // 5 kilos at 0 0 0 
        let mass1 = PointMass::new_gravity(
                Vec3::new(0_f32,0_f32,0_f32), // pos (at earth surface)
                Vec3::new(0_f32,0_f32,0_f32), // vel
                SimulationFloat::new(5_f64), // mass
                None, // model
            );
        // Three kilos 2 meters away from mass1
        let mass2 = PointMass::new_gravity(
                Vec3::new(2_f32,0_f32,0_f32), // pos (at earth surface)
                Vec3::new(1000_f32,1000_f32,1000_f32), // vel
                SimulationFloat::new(3_f64), // mass
                None, // model
            );
        let potential = mass1.calculate_gravitational_potential_energy_j(vec![ &mass2]);
        dbg!(potential.clone());
        assert!(potential < SimulationFloat::new(-5.005e-10_f32));
        assert!(potential > SimulationFloat::new(-5.006e-10_f32));
    }
    #[test]
    fn calc_potential_energy_different_orientation(){
        // 5 kilos at 1 1 0 
        let mass1 = PointMass::new_gravity(
                Vec3::new(1_f32,1_f32,0_f32), // pos (at earth surface)
                Vec3::new(0_f32,1_f32,0_f32), // vel
                SimulationFloat::new(5_f64), // mass
                None, // model
            );
        // Three kilos 2 meters away from mass1
        let mass2 = PointMass::new_gravity(
                Vec3::new(-1_f32,1_f32,0_f32), // pos (at earth surface)
                Vec3::new(1000_f32,1000_f32,1000_f32), // vel
                SimulationFloat::new(3_f64), // mass
                None, // model
            );
        let potential = mass1.calculate_gravitational_potential_energy_j(vec![ &mass2]);
        dbg!(potential.clone());
        assert!(potential < SimulationFloat::new(-5.005e-10_f32));
        assert!(potential > SimulationFloat::new(-5.006e-10_f32));
    }
    #[test]
    fn calc_kinetic_energy_test_zero(){
        // 5 kilos at 0 0 0 
        let mass = PointMass::new_gravity(
                Vec3::new(0_f32,0_f32,0_f32), // pos (at earth surface)
                Vec3::new(0_f32,0_f32,0_f32), // vel
                SimulationFloat::new(5_f64), // mass
                None, // model
            );
        assert_eq!(mass.calculate_kinetic_energy_j(), SimulationFloat::default());
    }
    #[test]
    fn calc_kinetic_energy_test_angles(){
        let speed = SimulationFloat::new(2_f32);
        for section in 1..10 {
            let mut angle1 = SimulationFloat::new(3.14159_f32 * 2_f32 / (section as f32)); 
            let mut angle2 = SimulationFloat::new(3.14159_f32 * 2_f32 / (section as f32)); 
            angle1.cos_mut();
            angle2.sin_mut();
            let x = speed.clone() * angle1;
            let y = speed.clone() * angle2;
            let mut mass = PointMass::new_gravity(
                    Vec3::new(0_f32,1_f32,0_f32), // pos 
                    Vec3::new(0_f32,0_f32,0_f32), // vel
                    SimulationFloat::new(5_f64), // mass
                    None, // model
                );
            mass.vel.x = x;
            mass.vel.y = y;
            let energy = mass.calculate_kinetic_energy_j();
            assert!(energy > SimulationFloat::new(9.99999999_f64));
            assert!(energy < SimulationFloat::new(10.0000001_f64));
        }
    }


}
