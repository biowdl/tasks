version 1.0
struct Readgroup {
    String id
    File R1
    String R1_md5
    File? R2
    String R2_md5
}

struct Library {
    String id
    Array[Readgroup]+ readgroups
}

struct Sample {
    String id
    Array[Library]+ libraries
}

task sampleConfigFileToStruct {
    input {
        File sampleConfigFile
        String outputJson = "output.json"
    }

    # Below command can convert any samplesheet with a nested dictionary
    # structure to a list of objects model.
    # It was specifically designed to run on both python2 and python3.
    # Only requirement is PyYAML.
    command {
        python <<CODE

        import yaml
        import json


        def nested_dicts_to_lists(dictionary):
            new_dict = dict()
            for key, value in dictionary.items():
                if type(value) == dict:
                    new_dict[key] = dict_to_item_list_with_id(value)
                else:
                    new_dict[key] = value
            return new_dict


        def dict_to_item_list_with_id(dictionary):
            items = []
            for sub_key, sub_dictionary in dictionary.items():
                item_dict = dict(id=sub_key, **nested_dicts_to_lists(sub_dictionary))
                items.append(item_dict)
            return items


        with open("~{sampleConfigFile}", "r") as samplesheet:
            samplesheet_dict = yaml.load(samplesheet)

        sample_struct = nested_dicts_to_lists(samplesheet_dict)

        with open("~{outputJson}", "w") as output_json:
            output_json.write(json.dumps(sample_struct))

        CODE
    }

    output {
        Map[String,Array[Sample]] map = read_json(outputJson)
        Array[Sample] samples = map["samples"]
     }
}