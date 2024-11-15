import pandas as pd
import numpy as np
import requests
import json
from typing import Iterable

from enums import ReferenceSpace


def clean_p4_section_images():
    # Read the P4 section images
    p4_section_images = pd.read_csv("data/P4_section_images.csv")

    # Extract the gene acronyms and section images
    gene_acronyms = p4_section_images.iloc[:, 0]
    section_images = p4_section_images.iloc[:, 1:]

    # Create a list to store the new rows
    new_rows = []

    # Iterate through the section images and gene acronyms
    for index, row in section_images.iterrows():
        gene_acronym = gene_acronyms[index]
        images = row.values

        # Iterate through the images
        for image in images:
            if pd.notna(image):
                # Create a new row and append it to the list
                new_rows.append({"Gene Acronym": gene_acronym, "Section Image": image})

        print(f"Progress: {index + 1}/{len(gene_acronyms)}")

    # Convert the list of new rows to a DataFrame
    new_data = pd.DataFrame(new_rows)

    # Save the new DataFrame to a CSV file
    new_data.to_csv("data/P4_section_images_cleaned.csv", index=False)


def get_p4_section_data_set_ids_set() -> set:
    # Read the cleaned P4 section images
    p4_section_images = pd.read_csv("data/P4_section_images_cleaned.csv")

    # Extract the section images
    section_images = p4_section_images.iloc[:, 1]

    # Create a set to store the section data set IDs
    section_data_set_ids_set = set()

    # Iterate through the section images
    for image in section_images:
        # Parse the image string as JSON

        # replace single quotes with double quotes
        image = image.replace("'", "\"")

        # replace the word None with null
        image = image.replace("None", "null")

        try:
            image = json.loads(image)

            # get the data_set_id
            data_set_id = image.get("data_set_id")
            if data_set_id is not None:
                section_data_set_ids_set.add(int(data_set_id))

        except json.JSONDecodeError as e:
            print("ERROR DECODING JSON: " + str(e))
            continue
        finally:
            continue

        # Add the section data set IDs to the set

    return section_data_set_ids_set


def get_p4_section_data_set_ids_list() -> list:
    # for each gene, get the section data set IDs
    p4_section_images = pd.read_csv("data/P4_section_images_cleaned.csv")

    section_data_set_ids_list = []

    # get unique gene acronyms
    gene_acronyms = p4_section_images["Gene Acronym"].unique()

    # also create a dataframe to store the results

    # call this endpoint for each gene
    #'https://api.brain-map.org/api/v2/data/SectionDataSet/query.json?criteria=[failed$eqfalse],products[id$eq3],genes[acronym$eq%27Abtb1%27]&include=genes,section_images,specimen(donor(age))'

    print(f"Found {len(gene_acronyms)} unique gene acronyms")
    progress = 0

    new_df = pd.DataFrame(columns=["gene_acronym", "section_data_set_id", "section_image"])

    gene_acronym_col = []
    section_data_set_id_col = []
    section_image_col = []

    while progress < len(gene_acronyms):
        gene_acronym = gene_acronyms[progress]
        try:
            url = f"https://api.brain-map.org/api/v2/data/SectionDataSet/query.json?criteria=[failed$eqfalse],products[id$eq3],genes[acronym$eq%27{gene_acronym}%27]&include=genes,section_images,specimen(donor(age))"

            req = requests.get(url)
            data = req.json()

            if data is None:
                continue

            # list of objects
            section_data_set_ids = data.get("msg")
            if section_data_set_ids is None:
                continue

            dataframes = []

            for section_data_set_id in section_data_set_ids:
                specimen = section_data_set_id.get("specimen")

                if specimen is None:
                    progress += 1
                    continue

                donor = specimen.get("donor")
                age = donor.get("age")

                if age is None:
                    progress += 1
                    continue

                days = age.get("days")

                if days is None:
                    progress += 1
                    continue

                if days == 4:
                    # this is a P4 section data set

                    # add the section data set ID to the list
                    # create a dataframe
                    # for each section image in the section images list, create a new row

                    gene_acronym_col.append(gene_acronym)
                    section_data_set_id_col.append(section_data_set_id.get("id"))
                    section_image_col.append(section_data_set_id.get("section_images"))

                    section_data_set_ids_list.append(section_data_set_id.get("section_images"))

                    """
                    pd = pd.DataFrame()
                    pd['gene_acronym'] = gene_acronym
                    pd['section_data_set_id'] = section_data_set_id.get("section_images").get("data_set_id")
                    pd['section_image'] = section_data_set_id.get("section_images").get("section_image")
                    section_data_set_ids_list.append(section_data_set_id.get("section_images"))
                    """

                    break

            # store the results in a dataframe

            section_data_set_ids_list.append(section_data_set_ids)

            progress += 1
            print(f"Progress: {progress}/{len(gene_acronyms)}")

        except Exception as e:
            print("ERROR: " + str(e))
            continue

    new_df['gene_acronym'] = gene_acronym_col
    new_df['section_data_set_id'] = section_data_set_id_col
    new_df['section_image'] = section_image_col

    # save the results to a CSV file
    new_df.to_csv("data/p4_section_data_set_ids_with_only_P4.csv", index=False)

    #

def get_section_images_from_allen_atlas(reference_space: ReferenceSpace = ReferenceSpace.P56, x: int = 0, y: int = 0,
                                        z: int = 0, section_data_set_ids: Iterable[int] = []) -> dict | None:
    # turn the section_data_set_ids into a string
    try:
        section_data_set_ids_str = ",".join(map(str, section_data_set_ids))

        url = f"http://api.brain-map.org/api/v2/reference_to_image/{reference_space.value}.json?x={x}&y={y}&z={z}&section_data_set_ids={section_data_set_ids_str}"

        req = requests.get(url)
        data = req.json()

        if data is None:
            return None

        msg = data.get("msg")
        if msg is None:
            return None

        return msg
    except Exception as e:
        print("ERROR: " + str(e))
        return None


def get_p4_xyz():
    # Read the P4 structure annotations
    p4_structure_annotations = pd.read_csv(
        "data/P4_structure_annotations.csv",
    )

    p4_structure_annotations.columns = [
        "A", "B", "C", "D", "X", "Y", "Z"
    ]

    # only keep X, Y, Z columns and remove the rest
    p4_structure_annotations = p4_structure_annotations[["X", "Y", "Z"]]

    return p4_structure_annotations


def get_p56_xyz():
    # Read the P4 structure annotations
    p56_structure_annotations = pd.read_csv(
        "data/P56_structure_annotations.csv",
    )

    p56_structure_annotations.columns = [
        "A", "B", "C", "X", "Y", "Z"
    ]

    # only keep X, Y, Z columns and remove the rest
    p56_structure_annotations = p56_structure_annotations[["X", "Y", "Z"]]

    return p56_structure_annotations


def main():
    """
    p4_annotations = get_p4_xyz()
    p56_annotations = get_p56_xyz()

    p4_section_datasets = get_p4_section_data_set_ids_set()

    print(p4_section_datasets)
    print(f"There are {len(p4_section_datasets)} unique section dataset IDs in P4 data")

    # get the first 10 section dataset IDs
    p4_section_datasets_sample = list(p4_section_datasets)[:1000]

    # at this point, there are 20,274 unique section dataset IDs in P4 data

    # for each voxel in p56, call the API to get the images from p4

    # create a dataframe to store the results
    results = pd.DataFrame(columns=["Section Dataset Id", "Section Image"])

    progress = 0

    try:
        for index, row in p56_annotations.iterrows():
            x = row["X"]
            y = row["Y"]
            z = row["Z"]

            section_data_set_ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            images = get_section_images_from_allen_atlas(reference_space=ReferenceSpace.P56, x=x, y=y, z=z,
                                                         section_data_set_ids=p4_section_datasets_sample)

            if images is not None:
                new_df = pd.DataFrame(images)
                results = pd.concat([results, new_df], ignore_index=True)
            else:
                print("No images found")

            progress += 1
            print(f"Progress: {progress}/{len(p56_annotations)}")
    except Exception as e:
        print("ERROR: " + str(e))

        # save the results to a CSV file
        results.to_csv("data/p4_images_from_p56_coords.csv", index=False)

    # save the results to a CSV file
    results.to_csv("data/p4_images_from_p56_coords.csv", index=False)

    pass
    """
    get_p4_section_data_set_ids_list()

if __name__ == "__main__":
    main()
