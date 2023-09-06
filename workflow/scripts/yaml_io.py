"""Yaml interactions"""
import yaml


def read_yaml(yaml_file_path: str) -> dict:
    with open(yaml_file_path, encoding="utf-8") as yaml_file:
        return yaml.load(yaml_file, Loader=yaml.FullLoader)


def write_yaml(yaml_file_path: str, dump_dict: dict) -> None:
    with open(yaml_file_path, "w", encoding="utf-8") as file:
        yaml.dump(dump_dict, file)
