import React from 'react';

export default function DownloadFile(props) {

  return (
    <p>
      <a href={props.file}>{props.name}</a>
    </p>
  )
}