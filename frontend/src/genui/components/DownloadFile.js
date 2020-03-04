import React from 'react';

export default function DownloadFile(props) {

  return (
    <a href={props.file}>{props.name}</a>
  )
}