function [ time_input ] = get_time_input( time_start, time_end, time_bounce, N, N_bounce )

time_bounce_floor = time_bounce - 0.5E-03;
time_bounce_ceil  = time_bounce + 0.5E-03;

N_prebounce  = int32( 0.25 * N * ( time_bounce_floor - time_start ) / ( time_end - time_start ) );
N_postbounce = N - N_prebounce - N_bounce + 1;

time_input_bounce     = linspace( time_bounce_floor, time_bounce_ceil, N_bounce );
time_input_prebounce  = linspace( time_start, time_bounce_floor, N_prebounce );
time_input_postbounce = linspace( time_bounce_ceil, time_end, N_postbounce );

time_input = unique( [ time_input_prebounce, time_bounce, time_input_bounce, time_input_postbounce ] );

end
