# Makefile para limpiar archivos del proyecto

# Variables
IMAGE_EXTENSIONS = *.png *.jpg *.jpeg *.gif *.bmp *.svg
DAT_FILES = *.dat

# Objetivo principal para limpiar todo
clean: clean-images clean-dat
	@echo "Limpieza completa terminada"

# Limpiar solo imágenes
clean-images:
	@echo "Eliminando archivos de imagen..."
	@for ext in $(IMAGE_EXTENSIONS); do \
		if ls $$ext 1> /dev/null 2>&1; then \
			rm -f $$ext; \
			echo "Eliminados: $$ext"; \
		fi; \
	done
	@echo "Limpieza de imágenes completada"

# Limpiar solo archivos .dat
clean-dat:
	@echo "Eliminando archivos .dat..."
	@if ls *.dat 1> /dev/null 2>&1; then \
		rm -f *.dat; \
		echo "Eliminados todos los archivos .dat"; \
	else \
		echo "No se encontraron archivos .dat"; \
	fi

# Limpiar archivos específicos de tu proyecto
clean-project:
	@echo "Eliminando archivos específicos del proyecto..."
	@rm -f analisis_errores.png
	@rm -f comparacion_soluciones_24.png
	@rm -f comparacion_soluciones_34.png
	@rm -f comparacion_soluciones.png
	@rm -f solucion_ejemplo1_analitica.dat
	@rm -f solucion_ejemplo1_numerica.dat
	@rm -f comparacion_numerica_24.png
	@rm -f solucion_numerica_34.png
	@echo "Archivos específicos del proyecto eliminados"

# Mostrar archivos que serían eliminados (modo dry-run)
preview:
	@echo "Archivos que serían eliminados:"
	@echo "=== IMÁGENES ==="
	@for ext in $(IMAGE_EXTENSIONS); do \
		if ls $$ext 1> /dev/null 2>&1; then \
			ls $$ext; \
		fi; \
	done 2>/dev/null || echo "No hay imágenes"
	@echo "=== ARCHIVOS .DAT ==="
	@ls *.dat 2>/dev/null || echo "No hay archivos .dat"

# Limpiar con confirmación
clean-confirm:
	@echo "¿Estás seguro de que quieres eliminar todas las imágenes y archivos .dat? [y/N]"
	@read -r response; \
	if [ "$$response" = "y" ] || [ "$$response" = "Y" ]; then \
		$(MAKE) clean; \
	else \
		echo "Operación cancelada"; \
	fi

# Ayuda
help:
	@echo "Makefile para limpiar archivos del proyecto"
	@echo ""
	@echo "Objetivos disponibles:"
	@echo "  clean          - Elimina todas las imágenes y archivos .dat"
	@echo "  clean images   - Elimina solo las imágenes"
	@echo "  clean dat      - Elimina solo los archivos .dat"
	@echo "  clean-project  - Elimina archivos específicos del proyecto"
	@echo "  preview        - Muestra qué archivos serían eliminados"
	@echo "  clean-confirm  - Limpia con confirmación del usuario"
	@echo "  help           - Muestra esta ayuda"

# Prevenir que make interprete los nombres como archivos
.PHONY: clean clean-images clean-dat clean-project preview clean-confirm help