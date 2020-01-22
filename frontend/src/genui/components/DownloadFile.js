import React from 'react';

export default function DownloadFile(props) {

  return (
    <p>
      <a href={props.file}a>{props.name}</a>
    </p>
  )
}