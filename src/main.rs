extern crate rand;

#[macro_use]
extern crate glium;
extern crate image;
extern crate time;

use std::io::Cursor;
use glium::{DisplayBuild, Surface};

use glium::glutin;

const MAZE_SIZE : i32 = 5;
const PLAYER_HEIGHT : i32 = 32;
const CELL_SIZE : i32 = 64;
const PROJPLANE_HEIGHT : i32 = 200;
const PROJPLANE_WIDTH : i32 = 320;
const PLAYER_FOV : f32 = 1.0472; // 60 degrees in radians.

// TODO: pack similar data into structs so it's easier to pass this stuff around...
fn draw_floor_column(pixel_x : i32,
                    pixel_y : i32,
                    bottom_y: i32,
                    wall_distance : f32,
                    wall_intersection_x : f32,
                    wall_intersection_y : f32,
                    distance_to_projplane : f32,
                    player_x : f32,
                    player_y : f32,
                    image_width : u32,
                    image_buffer : &Vec< Vec<(u8,u8,u8)> >,
                    dest_buffer : &mut Vec< Vec<(u8,u8,u8)> >) {
    // Floor Casting
    let mut current_pixel_y = pixel_y;
    let projection_plane_center_y = PROJPLANE_HEIGHT/2;

    while current_pixel_y >= bottom_y
    {
        let distance_from_center_to_proj_pixel = projection_plane_center_y - current_pixel_y;
        let staight_distance_from_player_to_floor_intersect = (PLAYER_HEIGHT as f32 / distance_from_center_to_proj_pixel as f32) * distance_to_projplane as f32;
        // Take the relative angle between the player_angle and the cast ray
        let distance_from_player_to_floor_intersect = staight_distance_from_player_to_floor_intersect; //* beta.cos(); TODO: get to work with proper distance.

        let weight = distance_from_player_to_floor_intersect / wall_distance;

        let floor_intersection_x = weight * wall_intersection_x + (1.0f32 - weight) * player_x;
        let floor_intersection_y = weight * wall_intersection_y + (1.0f32 - weight) * player_y;

        let floor_texture_offset_x = floor_intersection_x as u32 % CELL_SIZE as u32;
        let floor_texture_offset_y = floor_intersection_y as u32 % CELL_SIZE as u32;
        let tex_x = (floor_texture_offset_x as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;
        let tex_y = (floor_texture_offset_y as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

        dest_buffer[current_pixel_y as usize][pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        //     // Ceiling is the same, just mirrored.
        //     let ceiling_y = PROJPLANE_HEIGHT as usize - projplane_pixel_y as usize - 1usize; // Subtract one because it's zero indexed.
        //     dest_buffer[ceiling_y][projplane_pixel_x as usize] = image_buffer[tex_y as usize][tex_x as usize];

        current_pixel_y -= 1;
    }
}

fn draw_wall_column(column_height : i32,
                    column_bottom : i32,
                    tex_x : usize,
                    image_height : u32,
                    column : i32,
                    image_buffer : &Vec< Vec<(u8,u8,u8)> >,
                    dest_buffer : &mut Vec< Vec<(u8,u8,u8)> >) {
    for proj_y in column_bottom..column_height + column_bottom {
        // TODO: Look into this method of using ints to avoid floats.
        let d = proj_y * 256 - PROJPLANE_HEIGHT * 128 + column_height * 128;  //256 and 128 factors to avoid floats
        let tex_y = ((d * image_height as i32) / column_height) / 256;
        if proj_y >= 0 && proj_y < PROJPLANE_HEIGHT {
            dest_buffer[proj_y as usize][column as usize] = image_buffer[tex_y as usize % image_height as usize][tex_x];
        }
    }
}

fn main() {
    // building the display, ie. the main object
    let display = glutin::WindowBuilder::new()
        .with_vsync()
        .build_glium()
        .unwrap();

    // building a texture with "OpenGL" drawn on it
    let image = image::load(Cursor::new(&include_bytes!("floor_sample.png")[..]),
                            image::PNG).unwrap().to_rgba();
    let image_dimensions = image.dimensions();
    let (image_width, image_height) = image_dimensions;
    let mut image_raw = image.into_raw().clone();

    // Reversing the image data arranges the pixels correctly but reverses the rgba values of each so they must be reversed again when stored in the image buffer.
    image_raw.reverse();

    let mut image_buffer:Vec< Vec<(u8,u8,u8)> > = Vec::with_capacity(image_height as usize);

    for row_index in 0..image_height {
        let mut row = Vec::with_capacity(image_width as usize);

        for column_index in 0..image_width {

            let image_pixel = row_index as usize * (image_width * 4) as usize + (column_index * 4) as usize;
            //let image_pixel_a = image_raw[image_pixel];
            let image_pixel_b = image_raw[image_pixel+1];
            let image_pixel_g = image_raw[image_pixel+2];
            let image_pixel_r = image_raw[image_pixel+3];

            row.push((image_pixel_r,image_pixel_g, image_pixel_b));
        }

        image_buffer.push(row);
    }

    // building an empty texture
    let dest_texture = glium::Texture2d::empty_with_format(&display,
                                               glium::texture::UncompressedFloatFormat::U8U8U8U8,
                                               glium::texture::MipmapsOption::NoMipmap,
                                               PROJPLANE_WIDTH as u32, PROJPLANE_HEIGHT as u32).unwrap();

    dest_texture.as_surface().clear_color(0.0, 0.0, 0.0, 1.0);

    let mut dest_buffer:Vec< Vec<(u8,u8,u8)> > = Vec::with_capacity(PROJPLANE_HEIGHT as usize);

    for _ in 0..PROJPLANE_HEIGHT {
        let mut row = Vec::with_capacity(PROJPLANE_WIDTH as usize);

        // Init this pixel buffer to greeeeen!!
        for _ in 0..PROJPLANE_WIDTH {
            row.push((0u8,255u8,0u8));
        }

        dest_buffer.push(row);
    }

    // RAYCAST STUFFF!! PEW PEW
    let mut player_x = 160f32;
    let mut player_y = 160f32;
    let mut player_movement_forward = 0f32;
    let player_turn_speed = (360f32).to_radians();
    let player_move_speed = 100f32;
    let mut player_angle = (0f32).to_radians();

    let degree_between_columns = PLAYER_FOV / PROJPLANE_WIDTH as f32;

    let distance_to_projplane = (PROJPLANE_WIDTH as f32 / 2.0) / (PLAYER_FOV as f32 / 2.0).tan();

    // Contains the height of each floor tile.
    let maze = [64, 64, 64, 64, 64,// 0 - 64
                64, 0, 0, 32, 64, // 64 - 128
                64, 0, 0, 0, 56, // 128 - 196
                0, 64, 0, 0, 64, // 196 - 256
                0, 0, 32, 64, 0];


    //let mut last_time = time::precise_time_s() as f32;
    let target_frame_time = 1f32/60f32;
    let delta_time = target_frame_time;

    // Input
    let mut forward_button = false;
    let mut backward_buttom = false;
    let mut turn_left = false;
    let mut turn_right = false;

    'main: loop {
        { // Raycast rendering
            let player_cell_x = (player_x / CELL_SIZE as f32) as i32;
            let player_cell_y = (player_y / CELL_SIZE as f32) as i32;
            let current_height = maze[(player_cell_y * MAZE_SIZE + player_cell_x) as usize];

            // Start with the left side of the players view. Without the + 0.00001 there will be visual artifacts with rays that intersect walls completely perpendicular.
            let mut ray_angle : f32 = player_angle - (PLAYER_FOV / 2.0) + 0.00001;

            // Wrap
            if ray_angle < 0f32 {
                ray_angle = (360f32).to_radians() + ray_angle;
            }

            for column in 0..PROJPLANE_WIDTH {
                let ray_angle_degrees = ray_angle.to_degrees() as i32;

                // Initially the top of the last wall projected will be the bottom of the projection plane.
                let mut top_of_last_wall = 0i32;

                // horizontal wall intersection data for raycast.
                let mut horz_wall_intersection_y : f32;
                let horz_wall_y_offset : f32;
                let mut horz_wall_x_offset : f32;
                let mut horz_wall_intersection_x : f32;

                // vertical wall intersection data for raycast.
                let mut vert_wall_intersection_y : f32;
                let mut vert_wall_y_offset : f32;
                let vert_wall_x_offset : f32;
                let mut vert_wall_intersection_x : f32;

                horz_wall_x_offset = CELL_SIZE as f32 / ray_angle.tan();
                vert_wall_y_offset = CELL_SIZE as f32 * ray_angle.tan();

                // Facing up
                // Reminder coordinate system has y growing downward and x growing to the left
                if ray_angle >= (0f32).to_radians() && ray_angle < (180f32).to_radians() {
                    if vert_wall_y_offset < 0.0 {
                        vert_wall_y_offset *= -1f32;
                    }
                    // Point A will boarder the cell below
                    horz_wall_intersection_y = ((player_y as i32 / CELL_SIZE) * CELL_SIZE) as f32 + CELL_SIZE as f32;
                    horz_wall_y_offset = CELL_SIZE as f32;
                } else {
                    if vert_wall_y_offset > 0.0 {
                        vert_wall_y_offset *= -1f32;
                    }
                    //  Point A will boarder the cell above
                    horz_wall_intersection_y = ((player_y as i32 / CELL_SIZE) * CELL_SIZE) as f32;
                    horz_wall_y_offset = -CELL_SIZE as f32;
                }

                // TODO: Handle when the angle is zero/360 or 180 degrees (parallel to horizontal walls)
                let mut x_dist_to_intersection : f32 = std::f32::MAX;
                if ray_angle_degrees != 0 && ray_angle_degrees != 180 {
                    x_dist_to_intersection = (horz_wall_intersection_y - player_y) as f32 / ray_angle.tan();
                }

                horz_wall_intersection_x = x_dist_to_intersection + player_x as f32;

                // Facing right
                if ray_angle < (90f32).to_radians() || ray_angle >= (270f32).to_radians()
                {
                    if horz_wall_x_offset < 0f32 {
                        horz_wall_x_offset *= -1f32;
                    }
                    // Point A will boarder the cell to the right
                    vert_wall_intersection_x = ((player_x as i32 / CELL_SIZE) * CELL_SIZE) as f32 + CELL_SIZE as f32;
                    vert_wall_x_offset = CELL_SIZE as f32;
                }
                else
                {
                    if horz_wall_x_offset > 0f32 {
                        horz_wall_x_offset *= -1f32;
                    }
                    //  Point A will boarder the cell to the left
                    vert_wall_intersection_x = ((player_x as i32 / CELL_SIZE)  * CELL_SIZE) as f32; // TODO: try moving this -1 maybe it'll help fix the corner issue (def is having an effect).
                    vert_wall_x_offset = -CELL_SIZE as f32;
                }

                let y_dist_to_intersection : f32 = (vert_wall_intersection_x - player_x) as f32 * ray_angle.tan();
                vert_wall_intersection_y = y_dist_to_intersection + player_y as f32;
                // TODO: remove the rounded ray_angle_degrees
                if !(ray_angle >= (0f32).to_radians() && ray_angle < (180f32).to_radians()) {
                    horz_wall_intersection_y -= 1f32;
                }

                if !(ray_angle < (90f32).to_radians() || ray_angle >= (270f32).to_radians()) {
                    vert_wall_intersection_x -= 1f32;
                }

                // TODO: merge the vert and horz intersection tests into the same while loop! :)
                // TODO: do the raised/sunken checks with the heights in there and draw floors and walls as appropriate.

                /*let mut ray_casting = true;
                println!("starting");
                while ray_casting {
                        let mut horz_wall_height = std::f32::MAX;
                        let mut horz_wall_distance = std::f32::MAX;
                        let mut vert_wall_height = std::f32::MAX;
                        let mut vert_wall_distance = std::f32::MAX;
                        let mut horz_wall_intersected = false;
                        let mut vert_wall_intersected = false;
                        let mut horz_out_of_range = false;
                        let mut vert_out_of_range = false;

                        // Check if there has been a horizontal intersection with a wall.
                        // If the ray angle is angle is equal to 180 or 0, it will never intersect any horizontal wall.
                        if ray_angle_degrees != 180 && ray_angle_degrees != 0 {
                            let cell_y : i32 = (horz_wall_intersection_y / CELL_SIZE as f32) as i32;
                            let cell_x : i32 = (horz_wall_intersection_x / CELL_SIZE as f32) as i32;

                            horz_out_of_range = horz_wall_intersection_x < 0f32 || horz_wall_intersection_y < 0f32 || horz_wall_intersection_x >= (MAZE_SIZE * CELL_SIZE) as f32 || horz_wall_intersection_y >= (MAZE_SIZE * CELL_SIZE) as f32;

                            horz_wall_intersected = !horz_out_of_range && maze[(cell_y * MAZE_SIZE + cell_x) as usize] > 0;

                            if  horz_wall_intersected {
                                horz_wall_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                                let y_dist = (horz_wall_intersection_y - player_y).powf(2.0);
                                let x_dist = (horz_wall_intersection_x - player_x).powf(2.0);
                                horz_wall_distance = (x_dist + y_dist).sqrt();
                            }
                        }

                        if ray_angle.to_degrees() != 90f32 && ray_angle.to_degrees() != 270f32 {
                            // Check if there has been a vertical intersection with a wall.
                            let cell_y = (vert_wall_intersection_y / CELL_SIZE as f32) as i32;
                            let cell_x = (vert_wall_intersection_x / CELL_SIZE as f32) as i32;

                            vert_out_of_range = vert_wall_intersection_x < 0f32 || vert_wall_intersection_y < 0f32 || vert_wall_intersection_x >= (MAZE_SIZE * CELL_SIZE) as f32 || vert_wall_intersection_y >= (MAZE_SIZE * CELL_SIZE) as f32;
                            vert_wall_intersected = !vert_out_of_range && maze[(cell_y * MAZE_SIZE + cell_x) as usize] > 0;

                            if  vert_wall_intersected {
                                // Does the line of intersection belong to a cell that is a wall?
                                vert_wall_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                                let y_dist = (vert_wall_intersection_y - player_y).powf(2.0);
                                let x_dist = (vert_wall_intersection_x - player_x).powf(2.0);
                                vert_wall_distance = (x_dist + y_dist).sqrt();
                            }
                        }

                        if horz_wall_intersected || vert_wall_intersected
                        {
                            let mut height = horz_wall_height;
                            let mut distance = horz_wall_distance;
                            let mut wall_texture_offset = horz_wall_intersection_x as u32 % CELL_SIZE as u32;

                            if vert_wall_distance < horz_wall_distance
                            {
                                height = vert_wall_height;
                                distance = vert_wall_distance;
                                wall_texture_offset = vert_wall_intersection_y as u32 % CELL_SIZE as u32;
                            }

                            // Correct "fishbowl effect". Beta angle is the angle of the cast ray, relative to the player's viewing angle.
                            let beta = ray_angle - player_angle;
                            distance = distance * beta.cos();

                            let projected_slice_height = (height / distance * distance_to_projplane) as i32;
                            let projected_bottom_coord = PROJPLANE_HEIGHT/2 - (PLAYER_HEIGHT as f32 / (distance / distance_to_projplane)) as i32;
                            //let projected_bottom_coord = PROJPLANE_HEIGHT/2 - projected_slice_height/2;

                            let tex_x = (wall_texture_offset as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

                            draw_wall_column(projected_slice_height, projected_bottom_coord, tex_x, image_height, column, &image_buffer, &mut dest_buffer);
                        }

                        // If the rays checking horizontal and vertical intersections are both out of range, stop casting.
                        if horz_out_of_range || vert_out_of_range {
                            ray_casting = false;
                            continue;
                        }

                        // Extend the ray further using the intersection offsets, checking what will be the next intersection.
                        horz_wall_intersection_y += horz_wall_y_offset;
                        horz_wall_intersection_x += horz_wall_x_offset;

                        vert_wall_intersection_y += vert_wall_y_offset;
                        vert_wall_intersection_x += vert_wall_x_offset;
                }*/

                let mut dist_to_horz_intersection: f32 = std::f32::MAX;
                let mut hor_cell_height : f32 = 0.0;
                // If the ray angle is angle is equal to 180 or 0, it will never intersect any horizontal wall.
                if ray_angle_degrees != 180 && ray_angle_degrees != 0 {
                    // Will be used to index into the map array to check for walls.
                    let mut cell_y : i32 = (horz_wall_intersection_y / CELL_SIZE as f32) as i32;
                    let mut cell_x : i32 = (horz_wall_intersection_x / CELL_SIZE as f32) as i32;

                    let mut casting_ray = true;
                    while casting_ray {
                        // Check the horz_wall_intersections instead of the cell indexes for being < 0 since rounding (casting to i32) will sometimes round to 0 instead of a negative
                        // This will then display a horz wall at the 0 index outside the map.
                        if horz_wall_intersection_x < 0f32 || horz_wall_intersection_y < 0f32 || horz_wall_intersection_x >= (MAZE_SIZE * CELL_SIZE) as f32 || horz_wall_intersection_y >= (MAZE_SIZE * CELL_SIZE) as f32 {
                            dist_to_horz_intersection = std::f32::MAX;
                            casting_ray = false;
                        } else if maze[(cell_y * MAZE_SIZE + cell_x) as usize] > 0 {
                            hor_cell_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                            let y_dist = (horz_wall_intersection_y - player_y).powf(2.0);
                            let x_dist = (horz_wall_intersection_x - player_x).powf(2.0);
                            dist_to_horz_intersection = (x_dist + y_dist).sqrt();
                            casting_ray = false;
                        } else {
                            // Extend the ray further, checking what will be the next intersection.
                            horz_wall_intersection_y += horz_wall_y_offset;
                            horz_wall_intersection_x += horz_wall_x_offset;

                            cell_y = (horz_wall_intersection_y / CELL_SIZE as f32) as i32;
                            cell_x = (horz_wall_intersection_x / CELL_SIZE as f32) as i32;
                        }
                    }
                }

                let mut dist_to_vert_intersection: f32 = std::f32::MAX;
                let mut vert_cell_height : f32 = 0.0;

                if ray_angle.to_degrees() != 90f32 && ray_angle.to_degrees() != 270f32 {
                    // Will be used to index into the map array to check for walls.
                    let mut cell_y = (vert_wall_intersection_y / CELL_SIZE as f32) as i32;
                    let mut cell_x = (vert_wall_intersection_x / CELL_SIZE as f32) as i32;

                    let mut casting_ray = true;
                    while casting_ray {
                        // Check the vert_wall_intersections instead of the cell indexes for being < 0 since rounding (casting to i32) will sometimes round to 0 instead of a negative
                        // This will then display a vert wall at the 0 index outside the map.
                        if vert_wall_intersection_x < 0f32 || vert_wall_intersection_y < 0f32 || cell_y >= MAZE_SIZE || cell_x >= MAZE_SIZE {
                            dist_to_vert_intersection = std::f32::MAX;
                            casting_ray = false;
                        } else if maze[(cell_y * MAZE_SIZE + cell_x) as usize] > 0 {
                            // Does the line of intersection belong to a cell that is a wall?
                            vert_cell_height = maze[(cell_y * MAZE_SIZE + cell_x) as usize] as f32;
                            let y_dist = (vert_wall_intersection_y - player_y).powf(2.0);
                            let x_dist = (vert_wall_intersection_x - player_x).powf(2.0);
                            dist_to_vert_intersection = (x_dist + y_dist).sqrt();
                            casting_ray = false;
                        } else {
                            // Extend the ray further, checking what will be the next intersection.
                            vert_wall_intersection_y += vert_wall_y_offset;
                            vert_wall_intersection_x += vert_wall_x_offset;

                            cell_y = (vert_wall_intersection_y / CELL_SIZE as f32) as i32;
                            cell_x = (vert_wall_intersection_x / CELL_SIZE as f32) as i32;
                        }
                    }
                }

                let mut cell_height = vert_cell_height;
                let mut distance = dist_to_vert_intersection;
                let mut wall_texture_offset = vert_wall_intersection_y as u32 % CELL_SIZE as u32;
                let mut wall_intersection_x = vert_wall_intersection_x;
                let mut wall_intersection_y = vert_wall_intersection_y;

                if dist_to_horz_intersection < dist_to_vert_intersection {
                    cell_height = hor_cell_height;
                    distance = dist_to_horz_intersection;
                    wall_texture_offset = horz_wall_intersection_x as u32 % CELL_SIZE as u32;
                    wall_intersection_x = horz_wall_intersection_x;
                    wall_intersection_y = horz_wall_intersection_y;
                }

                // Correct "fishbowl effect". Beta angle is the angle of the cast ray, relative to the player's viewing angle.
                let beta = ray_angle - player_angle;
                distance = distance * beta.cos();

                let projected_slice_height = (cell_height / distance * distance_to_projplane) as i32;
                let projected_bottom_coord = PROJPLANE_HEIGHT/2 - (PLAYER_HEIGHT as f32 / (distance / distance_to_projplane)) as i32;
                //let projected_bottom_coord = PROJPLANE_HEIGHT/2 - projected_slice_height/2;

                let tex_x = (wall_texture_offset as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

                draw_wall_column(projected_slice_height, projected_bottom_coord, tex_x, image_height, column, &image_buffer, &mut dest_buffer);

                // Correct "fishbowl effect". Beta angle is the angle of the cast ray, relative to the player's viewing angle.
                /*let beta = ray_angle - player_angle;
                distance = distance * beta.cos();

                let projected_slice_height = (cell_height / distance * distance_to_projplane) as i32;
                let projected_bottom_coord = PROJPLANE_HEIGHT/2 - (PLAYER_HEIGHT as f32 / (distance / distance_to_projplane)) as i32;
                //let projected_bottom_coord = PROJPLANE_HEIGHT/2 - projected_slice_height/2;

                let tex_x = (wall_texture_offset as f32 * (image_width as f32 / CELL_SIZE as f32)) as usize;

                // Draw wall column
                for proj_y in projected_bottom_coord..projected_slice_height + projected_bottom_coord {
                    // TODO: Look into this method of using ints to avoid floats.
                    let d = proj_y * 256 - PROJPLANE_HEIGHT * 128 + projected_slice_height * 128;  //256 and 128 factors to avoid floats
                    let tex_y = ((d * image_height as i32) / projected_slice_height) / 256;
                    if proj_y >= 0 && proj_y < PROJPLANE_HEIGHT {
                        dest_buffer[proj_y as usize][column as usize] = image_buffer[tex_y as usize % image_height as usize][tex_x as usize];
                    }
                }*/

                // Only floor cast when a wall has been hit for now. 1000.0 is chosen at random pretty much.
                /*if distance < 1000.0
                {
                    // Floor Casting
                    let (projplane_pixel_x, projplane_pixel_y) = (column, projected_bottom_coord);

                    // TODO: Readd the ceilling drawing.
                    draw_floor_column(projplane_pixel_x,
                                    projplane_pixel_y,
                                    0,
                                    distance,
                                    wall_intersection_x,
                                    wall_intersection_y,
                                    distance_to_projplane,
                                    player_x,
                                    player_y,
                                    image_width,
                                    &image_buffer,
                                    &mut dest_buffer);
                } else {
                    println!("Oh no: {}, {}, angle: {}", vert_wall_intersection_x, vert_wall_intersection_y, ray_angle.to_degrees());
                }*/

                ray_angle += degree_between_columns;

                if ray_angle as i32 >= (360f32).to_radians() as i32 {
                    ray_angle = ray_angle - (360f32).to_radians();
                }
            }
        }

        // drawing a frame
        let target = display.draw();
        let the_texture = glium::Texture2d::new(&display, dest_buffer.clone()).unwrap();
        the_texture.as_surface().fill(&target, glium::uniforms::MagnifySamplerFilter::Nearest);
        target.finish().unwrap();

        // Clear Buffer
        dest_texture.as_surface().clear_color(0.0, 0.0, 0.0, 1.0);
        for row in 0..PROJPLANE_HEIGHT {
            for cell in 0..PROJPLANE_WIDTH {
                dest_buffer[row as usize][cell as usize] = (0u8, 255u8, 0u8);
            }
        }

        // TODO: Store the buttons which have been pressed in a map or array.
        for ev in display.poll_events() {
            match ev {
                glium::glutin::Event::Closed => return,
                glium::glutin::Event::KeyboardInput(element_state, _, virtual_key_code) => {
                    match element_state {
                        glium::glutin::ElementState::Pressed => {
                            let key = virtual_key_code.unwrap();
                            match key {
                                glium::glutin::VirtualKeyCode::W => {
                                    forward_button = true;
                                },
                                glium::glutin::VirtualKeyCode::A => {
                                    backward_buttom = true;
                                },
                                glium::glutin::VirtualKeyCode::S => {
                                    turn_left = true;
                                },
                                glium::glutin::VirtualKeyCode::D => {
                                    turn_right = true;
                                },
                                glium::glutin::VirtualKeyCode::Escape => {
                                    break 'main;
                                }
                                _ => ()
                            }
                        },
                        glium::glutin::ElementState::Released => {
                            let key = virtual_key_code.unwrap();
                            match key {
                                glium::glutin::VirtualKeyCode::W => {
                                    forward_button = false;
                                },
                                glium::glutin::VirtualKeyCode::A => {
                                    backward_buttom = false;
                                },
                                glium::glutin::VirtualKeyCode::S => {
                                    turn_left = false;
                                },
                                glium::glutin::VirtualKeyCode::D => {
                                    turn_right = false;
                                },
                                _ => ()
                            }
                        }
                    }
                },
                _ => ()
            }
        }

        if forward_button {
            player_movement_forward = player_move_speed;
        }
        if backward_buttom {
            player_angle -= player_turn_speed * delta_time;
        }
        if turn_left {
            player_movement_forward = -player_move_speed;
        }
        if turn_right {
            player_angle += player_turn_speed * delta_time;
        }

        player_x += player_movement_forward * delta_time * player_angle.cos();
        player_y += player_movement_forward * delta_time * player_angle.sin();
        player_movement_forward = 0f32;

        // Wrap
        if player_angle < 0f32 {
            player_angle = (360f32).to_radians() + player_angle;
        } else if player_angle >= (360f32).to_radians() {
            player_angle = player_angle - (360f32).to_radians();
        }

        std::thread::sleep(std::time::Duration::from_millis(16));
    }
}