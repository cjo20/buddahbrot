#ifdef WIN32 
#include <SDL.h>
#include <SDL_thread.h>
#else
#include <SDL2/SDL.h>
#include <SDL2/SDL_thread.h>
#endif
#include <stdio.h>
#include <cstdlib>
#include <complex>
#include "omp.h"
#include <random>
#include <chrono>
#include <thread>

#undef main

const int SCREEN_WIDTH = 1024;
const int SCREEN_HEIGHT = 1024;

const int RENDER_WIDTH = 4096;
const int RENDER_HEIGHT = 4096;

SDL_Window * window = NULL;
SDL_Surface * screen_surface = NULL;
SDL_Surface * hello_world = NULL;

bool quit = false;
bool waiting = false;
bool updated = true;

omp_lock_t writeLock;
omp_lock_t localLock;

typedef struct 
{
	double center_x;
	double center_y;
	double world_width;
	double world_height;
	int img_width;
	int img_height;
	unsigned char * pixels;
	std::mt19937 mrand;
} mandelbrot_data;

enum KeyPressSurfaces
{ 
	KEY_PRESS_SURFACE_DEFAULT, 
	KEY_PRESS_SURFACE_UP, 
	KEY_PRESS_SURFACE_DOWN, 
	KEY_PRESS_SURFACE_LEFT, 
	KEY_PRESS_SURFACE_RIGHT, 
	KEY_PRESS_SURFACE_TOTAL 
};

bool initWindow()
{
	bool success = true;

	if (SDL_Init(SDL_INIT_VIDEO) < 0)
	{
		printf("SDL could not initialize! SDL_ERROR: %s\n", SDL_GetError());
		success = false;
	}
	else
	{
		window = SDL_CreateWindow("Buddahbrot", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);

		if (window == NULL)
		{
			printf("Window could not be created! SDL_ERROR: %s\n", SDL_GetError());
			success = false;
		}
		else
		{
			screen_surface = SDL_GetWindowSurface(window);
		}
	}

	return success;
}

void closeWindow()
{
	SDL_DestroyWindow(window);
	SDL_Quit();
}

int mandelbrot(double x, double y, int max_iter, std::complex<double> * path)
{
	int i = 0;

	std::complex<double> c(x, y);
	std::complex<double> z = c;
	std::complex<double> old_z = 0;

	float p = sqrt(pow(x - 0.25, 2) + y*y);

	if (x < (p - 2*p*p + 0.25) && (pow(x+1, 2) + y*y) < 1/16.0)
	{
		return 0;
	}

	for (i = 0; i < max_iter && std::abs(z) < 2.0; ++i )
	{
		path[i] = z;
		if (z == old_z)
		{
			return 0;
		}

		if ((i & (i - 1)) == 0)
		{
			old_z = z;
		}

		z = z*z + c;
	}

	if (i == max_iter)
	{
		return 0;
	}
	else
	{
		return i;
	}

}

int buddahbrot(double epsilon, int max_iter, int image_width, unsigned * pixels, std::complex<double> * path, bool inv_y)
{
	int i = 0;

	int ix, iy;
	int bucket = 0;

	bool render_r = false;
	bool render_g = false;
	bool render_b = false;

#if 0
	if (max_iter < 50)
	{
		bucket = 0;
		return 0;
	}
	else if (max_iter < 500)
	{
		//return 0;
		render_r = true;
	}
	else if (max_iter < 4000)
	{
		render_b = true;
		render_r = true;
	}
	else
	{
		render_b = true;
	}
#else

	if (max_iter > 19 && max_iter < 5000)
	{
		render_b = true;
	}

	if (max_iter > 39 && max_iter < 10000)
	{
		render_g = true;
	}

	if (max_iter > 78 && max_iter < 20000)
	{
		render_r = true;
	}

#endif
	for (i = 0; i  < max_iter - 1; ++i)
	{
		ix = (real(path[i]) + 2.0) / epsilon;
		iy = (((inv_y? -1 : 1) *imag(path[i])) + 2.0) / epsilon;

		omp_set_lock(&localLock);
		if (render_r) pixels[3 * (iy * image_width + ix) + 0]++;
		if (render_g) pixels[3 * (iy * image_width + ix) + 1]++;
		if (render_b) pixels[3 * (iy * image_width + ix) + 2]++;
		omp_unset_lock(&localLock);
		
	}	
	return 1;
}


int buddahbrot_thread(void * ptr)
{
	mandelbrot_data * m_data = (mandelbrot_data *) ptr;	
	
	static unsigned pixels[RENDER_WIDTH * RENDER_HEIGHT * 3] = { 0 };
	
	const unsigned int max_iter = 20000;
	const double epsilon = 4.0 / RENDER_WIDTH;
	const unsigned int max_work = 2000000000;
	const unsigned int data_points = max_work / max_iter;

	std::complex<double> * path = new std::complex<double>[max_iter];
	std::uniform_real_distribution<double> dist(-2.0, 2.0);
	
	for (int r = 0; r < data_points; ++r)
	{
		double x, y;

		x = dist(m_data->mrand); //(mrand() * factor) - 2.0;
		y = dist(m_data->mrand); //(mrand() * factor) - 2.0;

		int iter = mandelbrot(x, y, max_iter, path);

		if (iter < max_iter && iter > 0)
		{
			buddahbrot(epsilon, iter, RENDER_WIDTH, pixels, path, false);
			buddahbrot(epsilon, iter, RENDER_WIDTH, pixels, path, true);
		}
	}

	omp_set_lock(&writeLock);
	{
		uint64_t total = 0;
		uint64_t count = 0;
		double average = 0.0;
		int max_val = 0;
		for (int x = 0; x < m_data->img_height * m_data->img_width * 3; ++x)
		{
			if (pixels[x] > 0)
			{
				total += pixels[x];
				count++;
			}
			if (pixels[x] > max_val)
			{
				max_val = pixels[x];
			}
		}

		average = total / (double) count;

		for (int x = 0; x < m_data->img_height * m_data->img_width * 3; ++x)
		{

			//unsigned char val = (unsigned char) ((255 * sqrt(pixels[x])) / sqrt(max_val / 4));
			//unsigned val = (unsigned char) ((255 * pixels[x]) / (max_val / 1));
			double val = ((255 * sqrt(pixels[x])) / sqrt(7 * average));

			if (val > 255)
			{
				val = 255;
			}

			m_data->pixels[x] = (unsigned char) val;
		}
	}
	omp_unset_lock(&writeLock);
	delete path;
	updated = true;

	return 0;	
}

int buddahbrot_runner(void * ptr)
{
	#pragma omp parallel shared(quit) num_threads(8)
	{
		int threadID = omp_get_thread_num();
		
		mandelbrot_data * m_data = (mandelbrot_data *) ptr;
		mandelbrot_data my_data;

		memcpy(&my_data, m_data, sizeof(mandelbrot_data));
		my_data.mrand = std::mt19937(omp_get_thread_num() * time(0));

		while (!quit)
		{
			buddahbrot_thread(&my_data);
		}
	}
	printf("Returning\n");
	return 0;
}

int main(int argc, char * argv[])
{
	unsigned char * pixels = new unsigned char [RENDER_WIDTH * RENDER_HEIGHT * 3];
	int i = 0;
	int writeCount = 0;
	SDL_Event e;
	mandelbrot_data m_data;

	omp_init_lock(&writeLock);
	omp_init_lock(&localLock);

	m_data.center_x = 0;
	m_data.center_y = 0;
	m_data.world_width = 1;
	m_data.world_height = 1;
	m_data.img_width = RENDER_WIDTH;
	m_data.img_height = RENDER_HEIGHT;
	m_data.pixels = pixels;

	SDL_Rect stretchRect;

	stretchRect.x = 0;
	stretchRect.y = 0;
	stretchRect.w = SCREEN_WIDTH;
	stretchRect.h = SCREEN_HEIGHT;

	initWindow();
	SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "2");
	SDL_Thread * threadID = SDL_CreateThread(buddahbrot_runner, "Mandelbrot Thread", &m_data);
	while (!quit)
	{
		while (SDL_PollEvent(&e) != 0)
		{
			if (e.type == SDL_QUIT)
			{
				quit = true;
			}
		}


		if (updated)
		{
			updated = false;
			hello_world = SDL_CreateRGBSurfaceFrom((void *) pixels, RENDER_WIDTH, RENDER_HEIGHT, 24, RENDER_WIDTH * 3, 0x0000FF, 0x00FF00, 0xFF0000, 0);
			SDL_BlitScaled(hello_world, NULL, screen_surface, &stretchRect);
			SDL_UpdateWindowSurface(window);
		}

		if (!(i % 200*10*10))
		{
			printf("Written %d times\n", ++writeCount);
			SDL_SaveBMP(hello_world, "Buddah.bmp");
		}	

		++i;
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}

	SDL_WaitThread(threadID, NULL);
	closeWindow();

	return 0;
}