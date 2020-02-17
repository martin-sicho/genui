import React from "react";
import { Input } from 'reactstrap';

function FileUpload(props) {
  const {field, form} = props;

  const handleChange = (e) => {
    e.preventDefault();
    const reader = new FileReader();
    const file  =  e.target.files[0];
    // const imgTag = document.getElementById("myimage");
    // imgTag.title = file.name;
    // reader.onload = function(event) {
    //   imgTag.src = event.target.result;
    // };
    reader.readAsDataURL(file);
    form.setFieldValue(field.name, file);
  };

  return (
    <Input
      name={field.name}
      id={field.name}
      type="file"
      onChange={handleChange}
    />
    // <div>
    //  {/*<input type={'file'} onChange={(o) => handleChange(o)} className={'form-control'}/>*/}
    //  {/*<img src={''} alt="" id={'myimage'}/>*/}
    // </div>
  );
}

export default FileUpload;