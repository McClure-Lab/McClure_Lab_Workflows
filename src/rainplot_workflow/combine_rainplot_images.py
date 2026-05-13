#!/usr/bin/env python3

import argparse
import os
from PIL import Image


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine rain plot PNGs with a genomic feature PNG."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to a text file containing rain plot PNG paths, one per line."
    )
    parser.add_argument(
        "--annotation_png",
        required=True,
        help="Path to the genomic feature PNG."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to write combined PNGs to."
    )
    parser.add_argument(
        "--delete_inputs",
        action="store_true",
        help="Delete the source rain plot and annotation PNGs after combining."
    )
    return parser.parse_args()


def read_manifest(manifest_path):
    with open(manifest_path, "r", encoding="utf-8") as handle:
        paths = [line.strip() for line in handle if line.strip()]
    return [path for path in paths if os.path.exists(path)]


def stack_images(top_path, bottom_path, output_path):
    with Image.open(top_path) as top_image, Image.open(bottom_path) as bottom_image:
        top = top_image.convert("RGB")
        bottom = bottom_image.convert("RGB")

        combined_width = max(top.width, bottom.width)
        combined_height = top.height + bottom.height

        canvas = Image.new("RGB", (combined_width, combined_height), "white")
        canvas.paste(top, ((combined_width - top.width) // 2, 0))
        canvas.paste(bottom, ((combined_width - bottom.width) // 2, top.height))
        canvas.save(output_path)


def resolve_annotation_path(rainplot_path, shared_annotation_path):
    annotation_name = os.path.basename(rainplot_path).replace(
        "rainplot_", "annotation_", 1
    )
    candidate = os.path.join(os.path.dirname(rainplot_path), annotation_name)
    if os.path.exists(candidate):
        return candidate
    return shared_annotation_path


def main():
    args = parse_args()
    rainplot_paths = read_manifest(args.manifest)

    if not rainplot_paths:
        raise FileNotFoundError(
            f"No rain plot PNGs were found from manifest: {args.manifest}"
        )

    if not os.path.exists(args.annotation_png):
        raise FileNotFoundError(
            f"Genomic feature PNG not found: {args.annotation_png}"
        )

    os.makedirs(args.output_dir, exist_ok=True)

    for rainplot_path in rainplot_paths:
        rainplot_name = os.path.basename(rainplot_path)
        combined_name = rainplot_name.replace("rainplot_", "combined_1", 1)
        output_path = os.path.join(args.output_dir, combined_name)
        annotation_path = resolve_annotation_path(rainplot_path, args.annotation_png)
        stack_images(rainplot_path, annotation_path, output_path)
        print(f"[INFO] Combined plot written to: {output_path}", flush=True)

        if args.delete_inputs and os.path.exists(rainplot_path):
            os.remove(rainplot_path)
        if (
            args.delete_inputs
            and annotation_path != args.annotation_png
            and os.path.exists(annotation_path)
        ):
            os.remove(annotation_path)

    if args.delete_inputs and os.path.exists(args.annotation_png):
        os.remove(args.annotation_png)


if __name__ == "__main__":
    main()
