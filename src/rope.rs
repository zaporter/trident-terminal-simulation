use std::collections::VecDeque;
use uuid::Uuid;
use crate::*;

#[derive(Debug, Clone)]
pub struct Rope{

    points : Vec<PointMass>,
    information_delay_s : SimulationFloat,
    state_queue : VecDeque<Vec<Vec3>>,
    total_simulated_time : SimulationFloat,
    force_uuid : Uuid,
    spring_constant : SimulationFloat,
    l_naught : SimulationFloat,
}
impl Rope {

    /*
     *
     * Beginning --- * --- * --- * --- End
     *
     *
     *
     */
    pub fn new(
            desired_length : SimulationFloat,
            num_springs : u32,
            density_kg_m : SimulationFloat,
            spring_constant : SimulationFloat, 
            speed_of_sound_m_s : SimulationFloat,
            beginning : PointMass, 
            end : PointMass,
            thickness_m : SimulationFloat,
            window:Option<&mut Window>,

        ) -> Rope {
        
        assert!(num_springs > 0);
        
        let difference = end.pos.clone() - beginning.pos.clone();
        let start_pos = beginning.pos.clone();
        let start_vel = beginning.vel.clone();
        let finish_vel = end.vel.clone();

        let mut points = Vec::new();

        let mut model : Option<SceneNode> = None;
        if let Some(window) = window {
            let mut model_unwrapped = window.add_sphere(thickness_m.to_f32());
            model_unwrapped.set_color(0.5,1.0,0.5);
            model = Some(model_unwrapped);
        }
        let rope_mass = desired_length.clone() * density_kg_m;
        for i in 1..num_springs {
            let percent_complete = SimulationFloat::new(i as f64/num_springs as f64);

            // This linearly interprets both the position and velocity from the beginning node to
            // the final node
            let point_pos = start_pos.clone() + (difference.clone() * percent_complete.clone());
            let point_vel = (start_vel.clone() * percent_complete.clone()) + (finish_vel.clone() * (SimulationFloat::new(1_f32)-percent_complete.clone()));

            let point = PointMass::new(
                point_pos,
                point_vel,
                rope_mass.clone() / SimulationFloat::new(num_springs), // mass,
                model.clone()
                );
            points.push(point);
               
        }
        points.push(beginning);
        points.push(end);
        let mut state_queue = VecDeque::new();
        let point_locations = points.iter().map(|p| p.pos.clone()).collect();
        state_queue.push_back(point_locations);
        Rope {
            points,
            information_delay_s : (desired_length.clone() / speed_of_sound_m_s) / SimulationFloat::new(num_springs),
            state_queue,
            total_simulated_time : SimulationFloat::default(),
            force_uuid : Uuid::new_v4(),
            spring_constant,
            l_naught : desired_length / SimulationFloat::new(num_springs),
        }

    }
    fn apply_external_forces(&mut self){
        // There should always be positions to use in the queue. 
        let positions : Vec<Vec3> = self.state_queue.pop_front().unwrap();
        let differences : Vec<Vec3>= positions.windows(2).map(|window| (window[0].clone() - window[1].clone())).collect();
        let forces : Vec<Vec3> = differences.iter().map(|d| d.clone().unit() *  (self.spring_constant.clone() * (d.modulus() - self.l_naught.clone()))).collect();
        


    }
}
impl Simulatable for Rope {
    fn tick(&mut self, dt_s : &SimulationFloat, gravity_applying_masses : &Vec<&PointMass>) {
        let minimum_end_time = self.total_simulated_time.clone() + dt_s.clone();
        // if self.total_simulated_time > minimum_end_time {
        //     for point in &mut self.points {
        //         point.tick(&self.information_delay_s, gravity_applying_masses);
        //     }
            
        // }
        while self.total_simulated_time < minimum_end_time {
            // update external forces with the state_queue
            for point in &mut self.points {
                point.tick(&self.information_delay_s, gravity_applying_masses);
            }
            self.total_simulated_time += self.information_delay_s.clone();
        }
    }
    fn calculate_momentum(&self) -> Vec3 {
        Vec3::default()
    }
    fn calculate_total_energy_j(&self, gravity_applying_masses : Vec<&PointMass>) -> SimulationFloat {
        SimulationFloat::default()

    }
    fn update_scene_node_position(&mut self){
        for point in &mut self.points {
            point.update_scene_node_position();
        }
    }
    fn draw_historical_positions(&self, window : &mut Window) {
        for point in &self.points{
            point.draw_historical_positions(window);
        }
    }

}
